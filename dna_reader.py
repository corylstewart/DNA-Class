import collections

class DnaFastQ(object):
    
    def __init__(self, filename):
        self.reads = None
        self.base_q = None
        self.filename = filename
        self.read_fastq_file()

    @staticmethod
    def encrypt_phred33(Q):
        return chr(int(Q)+33)

    @staticmethod
    def decrypt_phred33(letter):
        return ord(letter)-33

    def read_fastq_file(self):
        self.reads = []
        self.base_q = []
        with open(self.filename, 'r') as f:
            while True:
                if f.readline() == '':
                    break
                self.reads.append(f.readline().rstrip())
                f.readline()
                self.base_q.append([self.decrypt_phred33(x) for x in f.readline().rstrip()])
        return self.reads, self.base_q

    def find_gc_by_pos(self):
        gc = dict()
        totals = dict()
        for i in xrange(max([len(read) for read in self.reads])):
            gc[i] = 0.
            totals[i] = 0.
        for read in self.reads:
            for i in xrange(len(read)):
                totals[i] += 1
                if read[i] in ['G','C']:
                    gc[i] += 1
        return gc, totals

    def count_nucleotides(self):
        if not self.reads: return None
        else:
            totals = collections.Counter()
            for read in self.reads:
                totals.update(read)
        return totals


class DnaFa(object):

    def __init__(self, filename):
        self.genome = None
        self.filename = filename
        self.read_fa_file()

    def read_fa_file(self):
        with open(self.filename, 'r') as f:
            header = f.readline()
            self.genome = ''.join([snip.rstrip() for snip in f.readlines()])
        return self.genome

    def count_fa_freq(self):
        if self.fa_genome: return collections.Counter(self.fa_genome)
        else: return None