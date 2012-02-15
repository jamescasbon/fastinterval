"""
An simple interval class for DNA sequences

Typically, you will create a genome and then use that object to create
intervals.  The intervals have a sequence property that will look up the
actual sequence::

    >>> from fastinterval import Genome, Interval
    >>> test_genome = Genome('test/example.fa')
    >>> int1 = test_genome.interval(100, 150, chrom='1')
    >>> print int1
    1:100-150:
    >>> print int1.sequence
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA

fastinterval is using pyfasta to retrieve the sequence, so the access is mmapped.
It supports strandedness, which will be respected when accessing the sequence::

    >>> int2 = test_genome.interval(100, 150, chrom='1', strand=-1)
    >>> print int2.sequence
    TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC

The Interval class supports many interval operations::

    >>> int1 = test_genome.interval(100, 150, chrom='1')
    >>> int2 = test_genome.interval(125, 175, chrom='1')
    >>> int1.distance(int2)
    0
    >>> int1.span(int2)
    Interval(100, 175)
    >>> int1.overlaps(int2)
    True
    >>> int1.is_contiguous(int2)
    True
    >>> int1.contains(int2)
    False
    >>> int1.intersection(int2)
    Interval(125, 150)
    >>> int1.union(int2)
    Interval(100, 175)
    >>> Interval.merge([int1, int2, test_genome.interval(200,250, chrom='1')])
    [Interval(100, 175), Interval(200, 250)]

The Interval class is also based on bx python intervals.  So you can pass in
a value attritbue to point to an external object, and create interval trees and
so on.

    >>> from bx.intervals.intersection import IntervalTree
    >>> int3 = test_genome.interval(150, 200, chrom='1', value='foo')
    >>> tree = IntervalTree()
    >>> _ = map(tree.insert_interval, (int1, int2, int3))
    >>> tree.find(190, 195)
    [Interval(150, 200, value=foo)]

"""

VERSION = '0.0.1'

from pyfasta import Fasta
from bx.intervals import Interval as BaseInterval

def _convert_strand(strand):
    """ convert UCSC +/- to +1/-1"""
    if strand == '-': return -1
    if strand == '+': return 1
    return strand

class Interval(BaseInterval):
    """ A genomic interval """

    GENOME = None

    def __init__(self, start, stop, genome=None, **kws):
        self.genome = genome
        if 'strand' in kws:
            kws['strand'] = _convert_strand(kws['strand'])
        BaseInterval.__init__(self, start, stop, **kws)


    @property
    def sequence(self):
        """ Return the DNA from this intecal as a string """
        if not self.genome:
            raise Exception('Cannot retrieve sequence without a genome')

        return self.genome.sequence(dict(
            start = self.start,
            stop = self.end,
            chr = self.chrom,
            strand = self.strand
        ), one_based=False).upper()

    def __str__(self):
        return "%s:%s-%s:%s" % (self.chrom, self.start, self.end, self.strand if self.strand else '')

    def __len__(self):
        return self.end - self.start

    @classmethod
    def from_string(cls, loc, **kws):
        """ Create an interval from a chrx:start-end style string """
        # chr1:10858-10967:1
        toks = loc.split(':')
        chrom = toks[0]
        start, end = map(int, toks[1].split('-'))
        if len(toks) == 3:
            try:
                strand = int(toks[2])
            except ValueError:
                strand = toks[2]
                if strand == '+':
                    strand = 1
                elif strand == '-':
                    strand = -1
                elif strand == '':
                    strand = None
                else:
                    raise Exception('Unknown strand %s' % strand)
        else:
            strand = None
        return cls(start, end, chrom=chrom, strand=strand, **kws)

    def distance(self, other):
        """ return the distance between two intervals """
        if self.chrom != other.chrom:   return float('Inf')
        if self.start > other.end:      return self.start - other.end
        if self.end < other.start:      return other.start - self.end
        return 0

    def span(self, other, **kws):
        """ Return an interval spanning two intervals """
        if self.chrom != other.chrom:
            raise Exception('cannot get span over two chromosomes')
        return self.copy(
            start = min([self.start, other.start]),
            end = max([self.end, other.end]),
            **kws
        )

    def span_between(self, other, **kws):
        """Return an Inteval spanning the gap between two intervals """
        if self.chrom != other.chrom:
            raise Exception('cannot get span over two chromosomes')

        if self.overlaps(other):
            return None

        return self.copy(
            start = min([self.end, other.end]),
            end = max([self.start, other.start]),
            **kws
        )

    def overlaps(self, other):
        """ Return True if the intervals share at least one base """
        return self.distance(other) == 0 and not (self.start==other.end or self.end==other.start)

    def is_contiguous(self, other):
        """ Return True if the intervals are overlapping or contiguous """
        return self.distance(other) == 0

    def contains(self, other):
        """ Return true if one interval contains the other """
        return self.start <= other.start and self.end >= other.end

    def intersection(self, other):
        """ Return the interval containing the intersection of two intervals """
        if not self.overlaps(other):
            return None
        return Interval(
            max(self.start, other.start),
            min(self.end, other.end),
            chrom = self.chrom
        )

    def union(self, other, merge_contiguous=False):
        """ Return an interval containing the interval of two overlapping interval"""
        test = Interval.is_contiguous if merge_contiguous else Interval.overlaps
        if not test(self, other):
            raise Exception('cannot get union of non overlapping intervals')
        return self.span(other)

    def copy(self, **kws):
        """ Copy this interval, and optionally provide a dict of new attrs """
        if 'start' in kws:
            start = kws.pop('start')
        else:
            start = self.start
        if 'end' in kws:
            end = kws.pop('end')
        else :
            end = self.end

        template = dict(chrom=self.chrom, strand=self.strand, value=self.value)
        template.update(kws)

        return Interval(start, end, **template)


    def __sub__(self, other):
        """ Subtract an interval and return a list of intervals """
        if not self.overlaps(other):
            return [self]

        if other.contains(self):
            return []

        elif self.contains(other) and not (self.start == other.start or self.end == other.end):
            left = self.copy()
            left.end = other.start
            right = self.copy()
            right.start = other.end
            return [left, right]

        elif self.start >= other.start:
            right = self.copy()
            right.start = other.end
            return [right]

        elif self.end <= other.end:
            left = self.copy()
            left.end = other.start
            return [left]

        raise Exception('unhandled subtraction')



    @classmethod
    def merge(cls, intervals, merge_contiguous=False, **kwargs):
        """ merge a list of intervals and return a list of intervals

        By default, the intervals must be overlapping to be merged.  If you
        want to merge contiguous intervals, set merge_contiguous to True.
        """

        is_overlapping = cls.is_contiguous if merge_contiguous else cls.overlaps
        intervals = sorted(intervals, key=lambda x: (x.chrom, x.end))

        if len(intervals) > 1:
            done, todo = [intervals[0].copy(**kwargs)], intervals[1:]
            while todo:
                item = todo.pop(0).copy(**kwargs)

                while done and is_overlapping(item, done[-1]):
                    item = done.pop().union(item, merge_contiguous=merge_contiguous)
                done.append(item)
        else:
            done = intervals
        return done

    @classmethod
    def coverage(cls, intervals):
        if not intervals: return []
        chrom = intervals[0].chrom
        starts = [x.start for x in intervals]
        ends = [x.end for x in intervals]
        breaks = sorted(set(starts + ends))

        return [
            cls(start, end, chrom=chrom,
                value=1 + len([x for x in starts if x<= start]) - len([x for x in ends if x <= end])
            )
            for start, end in zip(breaks, breaks[1:])
        ]

    def add_border(self, size=0, upstream=0, downstream=0):
        """ return interval with some bases added to each end """
        if not (size or upstream or downstream):
            return self.copy()

        if size and (upstream or downstream):
            raise Exception('please either size or upstream/downstream')

        if size:
            return self.copy(start=self.start-size, end = self.end+size)

        if not self.strand:
            raise Exception('Cannot add upstrea/downstream to strandless interval')

        if self.strand == 1:
            return self.copy(start=self.start - upstream, end = self.end + downstream)

        else:
            return self.copy(start=self.start - downstream, end = self.end + upstream)

    def truncate(self, size):
        """ truncate this interval to size, respecting the orientation """
        if self.strand is None:
            raise Exception('cannot truncate unstranded interval')
        elif self.strand > 0:
            return self.copy(end=min(self.start+size, self.end))
        elif self.strand < 0:
            return self.copy(start=max(self.end-size, self.start))


class Genome(object):
    """ A convienience for creating intervals on the same genome """

    def __init__(self, fname, *args, **kws):
        """ Create a genome using a Fasta file. Other args passed to pyfasta """
        self.fasta = Fasta(fname, *args, **kws)

    def interval(self, start, end, **kws):
        """ return an interval on this genome """
        return Interval(start, end, genome=self.fasta, **kws)



class MinimalSpanningSet(object):
    """ Create a minimal spanning set for targets from a set of candidates """

    def score_candidate(self, candidate):
        return sum([
            len(candidate.intersection(x))
            for x in self.remaining_targets
            if candidate.overlaps(x)
        ])

    def coverage(self, choices):
        """ work out the covered based for a set of choices"""
        covered = Interval.merge(
            [
                choice.intersection(x)
                for x in self.targets
                for choice in choices
                if choice.overlaps(x)
            ]
        )
        return sum(map(len, covered))


    def __init__(self, targets, candidates, score_function=None, sort_key=None):
        self.targets = targets
        self.remaining_targets = list(targets)
        self.candidates = candidates
        self.chosen = []
        self.sort_key = sort_key
        if score_function is None:
            self.score_function = MinimalSpanningSet.score_candidate
        self._find_set()


    def _find_set(self):
        """ main loop """

        while True:
            scores = {}

            # work out the score by summing the length of overlaps with the target
            for candidate in self.candidates:
                scores[candidate] = self.score_candidate(candidate)

            # choose the best candidates
            rankings = sorted(self.candidates, key=scores.get, reverse=True)
            if rankings == []:
                break

            best = rankings[0]

            # break if no improvement is possible
            if scores[best] == 0:
                break

            if self.sort_key:
                equiv = [x for x in rankings if scores.get(x)==scores.get(best)]
                equiv = sorted(equiv, key=self.sort_key)
                best = equiv[-1]


            # update the data
            self.chosen.append(best)
            self.candidates.remove(best)
            self._update_targets(best)

            # break if no targets or candidates left
            if len(self.targets) == 0 or len(self.candidates) == 0:
                break

        self._remove_redundant()

    def _update_targets(self, choice):
        """ remove chosen interval from targets """
        new_targets = []
        for target in self.remaining_targets:
            new_targets.extend(target - choice)
        self.remaining_targets = new_targets

    def _remove_redundant(self):
        """ drop any candidates that are completely covered by the rest"""

        total_coverage = self.coverage(self.chosen)

        for c in list(self.chosen):
            # if the coverage without a choice is the same, remove it
            test = list(self.chosen)
            test.remove(c)
            if self.coverage(test) == total_coverage:
                self.chosen.remove(c)



