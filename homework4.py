from dna_reader import DnaFastQ
import itertools
import pickle
import time
from bisect import *


#start bisect add ons
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_gt(a, x):
    'Find leftmost value greater than x'
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError
#end bisect add ons

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def join_reads(a, b, o):
    return a + b[o:]

def preprocess(reads, min_overlap):
    d = dict()
    for i in range(len(reads)):
        for j in range(i+1, len(reads)):
            a = overlap(reads[i], reads[j], 0)
            if a >= min_overlap:
                d[(reads[i], reads[j])] = a
            a = overlap(reads[j], reads[i], 0)
            if a >= min_overlap:
                d[(reads[j], reads[i])] = a
    return d


def save_my_pre(filename, p_p):
    with open(filename, 'w+') as f:
        f.write(pickle.dumps(p_p))

def load_pre_proc_file(filename):
    with open(filename, 'r') as f:
        pre_pickle = f.read()
    return pickle.loads(pre_pickle)

def make_ordered_overlap(pre_proc):
    overlaps = list()
    for key, value in pre_proc.items():
        overlaps.append((-value, key))
    overlaps.sort()
    return overlaps

def scs(reads, o_o, min_ol=1):
    for i in range(len(o_o)):
        if o_o[i][0] > -min_ol:
            break
    old_oo = o_o[:]
    o_o = o_o[:i]
    reads.sort()
    read_dict = dict()
    for read in reads:
        read_dict[read] = 1
    read_len = len(reads[0])
    remove_matching_pairs(reads, o_o)
    while o_o and len(reads) > 1:
        process_next_pair(reads, read_dict, o_o, min_ol)
    return reads
            
def remove_matching_pairs(reads, o_o):
    if o_o[0][1][0] != o_o[0][1][1]:
        return
    c_pair = o_o.pop(0)
    while c_pair[1][0] == c_pair[1][1]:
        ind = index(reads, c_pair[1][0])
        while reads[ind] == reads[ind+1]:
            reads.pop(ind+1)
        c_pair = o_o.pop(0)
    insort(o_o, c_pair)

def process_next_pair(reads, read_dict, o_o, min_ol):
    remove_matching_pairs(reads, o_o)
    c_pair = o_o.pop(0)
    o_l = -c_pair[0]
    one = c_pair[1][0]
    two = c_pair[1][1]

    if one in read_dict: one_in = True
    else: one_in = False

    if two in read_dict: two_in = True
    else: two_in = False

    if one_in and two_in:
        add_new_read(reads, read_dict, o_o, one, two, o_l, min_ol)

def add_new_read(reads, read_dict, o_o, one, two, o_l, min_ol):
    read_dict.pop(one, None)
    read_dict.pop(two, None)
    reads.pop(index(reads, one))
    reads.pop(index(reads, two))
    new_read = join_reads(one, two, o_l)
    for read in reads:
        suf = overlap(read, new_read, 0)
        if suf >= min_ol:
            insort(o_o, (-suf, (read, new_read)))
        pre = overlap(new_read, read, 0)
        if pre >= min_ol:
            insort(o_o, (-pre, (new_read, read)))
    read_dict[new_read] = 1
    insort(reads, new_read)

def count_letter(s, l):
    count = 0
    for letter in s:
        if letter == l:
            count += 1
    return count

print 'starting'

#p_p = preprocess(reads, 1) #preprocess file filter out 0 overlaps
#ordered_overlap = make_ordered_overlap(p_p) #sort the dict
#save_my_pre(filename, ordered_overlap) #save preprocess for later use
start = time.time()
filename = 'data/pre_proc_reads.pickle'
reads = DnaFastQ('data/week4.fq').reads
ordered_overlap = load_pre_proc_file(filename)
print 'load data', time.time()-start

start = time.time()
z = scs(reads[:], ordered_overlap[:], 20)
print len(z[0]), time.time()-start, len(z)

print 'done'