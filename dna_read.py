from dna_reader import DnaFastQ

class DnaRead(object):
    def __init__(self, read):
        self.read = read

    def make_k_mers(self, k):
        k_mers = set()
        for i in range(len(self.read)-k+1):
            k_mers.add(self.read[i:i+k])
        k_mers = list(k_mers)
        #k_mers.sort()
        return k_mers

class OverlapFinder(object):

    def __init__(self, filename, k):
        self.k = k
        self.reads = self.generate_reads(filename)
        self.read_directory = self.make_read_directory()
        self.overlaps = 0
        self.overlap_graph_out, self.overlap_graph_in, \
            self.no_overlap = self.make_overlap_graph()

    def generate_reads(self, filename):
        fast_q = DnaFastQ(filename)
        reads = list()
        for read in fast_q.reads:
            reads.append(DnaRead(read))
        return reads

    def make_read_directory(self):
        d = dict()
        for read in self.reads:
            k_mers = read.make_k_mers(self.k)
            for k_mer in k_mers:
                if k_mer not in d:
                    d[k_mer] = list()
                d[k_mer].append(read)
        return d

    def overlap(self, read_one, read_two):
        a = read_one.read
        b = read_two.read
        start = 0  # start all the way at the left
        while True:
            start = a.find(b[:self.k], start)  # look for b's prefix in a
            if start == -1:  # no more occurrences to right
                return 0
            # found occurrence; check for full suffix/prefix match
            if b.startswith(a[start:]):
                return len(a)-start
            start += 1  # move just past previous match

    def make_overlap_graph(self):
        ol_g = dict()
        in_g = dict()
        no_overlap = list()
        for read_one in self.reads:
            best_overlap = 0
            best_read = None
            suf = read_one.read[-self.k:]
            for read_two in self.read_directory[suf]:
                if read_one != read_two:
                    overlap = self.overlap(read_one, read_two)
                    if overlap > 0:
                        self.overlaps += 1
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best_read = read_two
            if best_read:
                if read_one not in ol_g:
                    ol_g[read_one] = dict()
                if best_read not in in_g:
                    in_g[best_read] = dict()
                ol_g[read_one][best_read] = 1
                in_g[best_read][read_one] = 1
            else:
                no_overlap.append(read_one)
        return ol_g, in_g, no_overlap