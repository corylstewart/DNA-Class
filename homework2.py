from kmer_index import *
from dna_reader import DnaFa
from pigeon_hole import PigeonHole
from word_search import WordSearch
import bm_preproc 

print 'starting'
t = DnaFa('data/chr1.GRCh38.excerpt.fasta').genome

print 'Solving Q1'
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
one = len(t) - len(p) + 1
print 'Answer Q1', one
print ''

print 'Solving Q2'
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
two = WordSearch().naive_find_matches_count(p,t)
print 'Answer Q2', two
print ''

print 'Solving Q3'
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
p_bm = bm_preproc.BoyerMoore(p)
three = bm_preproc.boyer_moore_with_counter(p, p_bm, t)
print 'Answer Q3', three
print ''

print 'Solving Q4'
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
four = len(PigeonHole().approximate_matchs(p, t, 2))
print 'Answer Q4', four
print ''

print 'Solving Q5'
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
five = PigeonHole().approximate_matchs_with_count(p, t, 2)
print 'Answer Q5', five
print ''


print 'Solving Q6'
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
six = PigeonHole().approximate_matchs_using_subseq(p, t, 2, 3)
print 'Answer Q6', six