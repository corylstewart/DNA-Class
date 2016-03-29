import edit_distance as ED
from dna_reader import DnaFa
from dna_read import OverlapFinder

print 'starting'

def my_dist(p, t):
    len_p = len(p)
    len_t = len(t)
    best = len_p
    for i in range(len_t-len_p):
        if i%10000 == 0:
            print i, ' of ', len_t, ' me ', best
        for j in range(max(0, len_p-best),min(len_t, i + len_p + best)):
            #print max(0, len_p-best), min(len_t, i + len_p + best), len_p, best
            #raw_input()
            new_t = t[i:j]
            me = ED.EditDistance(p, new_t).make_edit_distance()
            if me < best:
                best = me
        if best == 0:
            return 0
    return best


t = DnaFa('data/chr1.GRCh38.excerpt.fasta').genome



print 'Start Q1'
p = 'GCTGATCGATCGTACG'
#one = ED.EditDistance(p, t).calculate_min_within_large_string() # 3
one = 3
print 'Solution Q1', one
print ''

print 'Start Q2'
p = 'GATTTACCAGATTGAG'
#two = ED.EditDistance(p, t).calculate_min_within_large_string() # 2
two = 2
print 'Solution Q2', two
print ''

print 'Start Q3'
of = OverlapFinder('data/week_3.fastq', 30)
print 'Solution Q3', of.overlaps
print ''

print 'Start Q4'
print 'Solution Q4', len(of.overlap_graph_out)