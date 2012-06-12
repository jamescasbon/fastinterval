from fastinterval import Interval, Genome, MinimalSpanningSet
import fastinterval
import pyfasta
import doctest


suite = doctest.DocTestSuite(fastinterval)

def test_Interval_from_string():
    a = Interval.from_string('chr1:10858-10967:1')
    assert a.chrom == 'chr1'
    assert a.start == 10858
    assert a.end == 10967
    assert a.strand == 1

    a = Interval.from_string('chr1:10858-10967:-1')
    assert a.chrom == 'chr1'
    assert a.start == 10858
    assert a.end == 10967
    assert a.strand == -1

def test_Interval_sequence():
    genome = pyfasta.Fasta('test/example.fa')
    l1 = Interval.from_string('1:858-967:1', genome=genome)
    l2 = Interval.from_string('1:858-967:-1', genome=genome)
    print l1.sequence
    print l2.sequence
    assert l1.sequence != l2.sequence


def test_Genome():
    genome = Genome('test/example.fa')
    i1 = genome.interval(10, 20, chrom='1')
    print 'genome is', i1.genome
    assert i1.genome == genome.fasta
    assert i1.sequence


def test_Interval_distance():
    l1 = Interval.from_string('chr1:10858-10967:1')
    l2 = Interval.from_string('chr1:10858-10967:-1')
    assert l1.distance(l2) == 0
    l3 = Interval.from_string('chr1:10968-10977:-')
    print l1.distance(l3)
    assert l1.distance(l3) == 1
    assert l3.distance(l1) == 1

def test_Interval_span():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:20858-30000:-1')
    l3 = Interval.from_string('chr2:20858-30000:-1')
    ex = l1.span(l2)
    print ex
    assert ex.start == 10000
    assert ex.end == 30000
    ex = l2.span(l1)
    assert ex.start == 10000
    assert ex.end == 30000
    try:
        l1.span(l3)
        assert False
    except:
        pass

def test_Interval_overlaps():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-21000:-1')
    l3 = Interval.from_string('chr1:20858-30000:-1')
    l4 = Interval.from_string('chr1:30000-30001:-1')
    assert l1.overlaps(l2)
    assert l2.overlaps(l3)
    assert not l1.overlaps(l3)
    assert not l3.overlaps(l4)

def test_Interval_is_contiguous():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-21000:-1')
    l3 = Interval.from_string('chr1:20858-30000:-1')
    l4 = Interval.from_string('chr1:30000-30001:-1')
    print l1.distance(l2) == 0
    assert l1.is_contiguous(l2)
    assert l2.is_contiguous(l3)
    assert not l1.is_contiguous(l3)
    assert  l3.is_contiguous(l4)

def test_Interval_contains():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-10964:-1')
    l3 = Interval.from_string('chr1:20858-30001:-1')
    l4 = Interval.from_string('chr1:30000-30003:-1')
    assert l2 in l1
    assert not l1 in l2
    assert not l3 in l4

def test_Interval_intersection():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-10964:-1')
    i = l1.intersection(l2)
    assert i.start == l2.start
    assert i.end == l2.end
    assert i.chrom == l2.chrom

    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-11000:-1')
    i = l1.intersection(l2)
    assert i.start == l2.start
    assert i.end == l1.end
    assert i.chrom == l2.chrom

    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr2:10858-10964:-1')
    i = l1.intersection(l2)
    assert i is None

def test_Interval_union():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-10964:-1')
    i = l1.union(l2)
    assert i.start == l1.start
    assert i.end == l1.end
    assert i.chrom == l2.chrom

    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-11000:-1')
    i = l1.union(l2)
    assert i.start == l1.start
    assert i.end == l2.end
    assert i.chrom == l2.chrom

    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr2:10858-10964:-1')
    i = l1.intersection(l2)
    assert i is None

def test_Interval_sub():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-10964:-1')
    l3 = Interval.from_string('chr2:10858-10964:-1')
    i = l1 - l2
    print i
    assert len(i)==2
    l = i[0]
    r = i[1]

    assert l.start == 10000
    assert l.end == 10858
    assert r.start == 10964
    assert r.end == 10967

    assert len(l1 - l3) == 1

def test_Interval_merge():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = Interval.from_string('chr1:10858-12964:-1')
    l3 = Interval.from_string('chr1:10858-10964:-1')

    merged = Interval.merge([l1,l2,l3])
    assert len(merged) == 1
    m = merged[0]
    assert m.start == 10000
    assert m.end == 12964

def test_find_minimal_spanning_set():

    targets = [
        Interval.from_string('chr1:10000-11000:+'),
        Interval.from_string('chr1:12000-13000:+')
    ]

    candidates = [
        Interval.from_string('chr1:10000-10500:+'),
        Interval.from_string('chr1:10100-10600:+'),
        Interval.from_string('chr1:10300-11700:+'),
        Interval.from_string('chr1:10500-11000:+'),
        Interval.from_string('chr1:10700-11100:+'),
    ]

    reads = MinimalSpanningSet(targets, candidates)
    n_reads = len(reads.chosen)
    print 'chose', reads.chosen
    print 'n_reads', n_reads
    assert n_reads == 2

    targets = [
        Interval.from_string('chr1:12000-13000:+')
    ]

    candidates = [
        Interval.from_string('chr1:10000-10500:+'),
        Interval.from_string('chr1:10100-10600:+'),
        Interval.from_string('chr1:10300-11700:+'),
        Interval.from_string('chr1:10500-11000:+'),
        Interval.from_string('chr1:10700-11100:+'),
    ]


    # this example is constructed into fooling the greedy algorithm
    # to choose one too many
    reads = MinimalSpanningSet(targets, candidates)
    assert len(reads.chosen) == 0
    targets = [
        Interval.from_string('chr1:12000-13000:+')
    ]

    candidates = [
        Interval.from_string('chr1:12000-12500:+'),
        Interval.from_string('chr1:12500-13000:+'),
        Interval.from_string('chr1:12100-12900:+'),
    ]

    reads = MinimalSpanningSet(targets, candidates)
    assert len(reads.chosen) == 2

def test_add_border():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = l1.add_border(upstream=50, downstream=100)
    assert l2.start == 10000 - 50
    assert l2.end == 10967 + 100

def test_truncate():
    l1 = Interval.from_string('chr1:10000-10967:1')
    l2 = l1.truncate(100)
    assert l2.start == 10000
    assert l2.end == 10100

    l1 = Interval.from_string('chr1:10000-10967:-1')
    l2 = l1.truncate(100)
    assert l2.start == 10867
    assert l2.end == 10967



