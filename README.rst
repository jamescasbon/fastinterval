
An simple interval class for DNA sequences from FASTA files that provides
fast access to sequences and methods for interval logic on those sequences.

An online version of the documentation should be available at 
http://fastinterval.readthedocs.org/

Usually you will create a `Genome` and then use that object to create
intervals.  The intervals have a sequence property that will look up the
actual sequence::

    >>> from fastinterval import Genome, Interval
    >>> test_genome = Genome('test/example.fa')
    >>> int1 = test_genome.interval(100, 150, chrom='1')
    >>> print int1
    1:100-150:
    >>> print int1.sequence
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA

fastinterval uses pyfasta to retrieve the sequence, so the access is mmapped
(i.e fast).  It supports strandedness, which will be respected when accessing
the sequence::

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
    >>> int1 in int2
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

