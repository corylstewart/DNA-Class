from kmer_index import * 

class PigeonHole(object):

    @staticmethod
    def approximate_matchs(p, t, n):
        segment_length = int(round(len(p)/(n+1)))
        t_index = Index(t,segment_length)
        all_matches = set()
        for i in range(n+1):
            start = i*segment_length
            end = min(start+segment_length, len(p))
            matches = t_index.query(p[start:end])
            for m in matches:
                if m < start or m-start+len(p) > len(t):
                    continue
                mismatches = 0
                for j in range(start):
                    if p[j] != t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break
                    if mismatches > n:
                        break
                for j in range(end, len(p)):
                    if p[j] != t[m-start+j]:
                        mismatches += 1
                        if mismatches > n:
                            break
                if mismatches <= n:
                    all_matches.add(m-start)
        return all_matches

    @staticmethod
    def approximate_matchs_with_count(p, t, n):
        segment_length = int(round(len(p)/(n+1)))
        t_index = Index(t,segment_length)
        all_matches = set()
        hits = 0
        for i in range(n+1):
            start = i*segment_length
            end = min(start+segment_length, len(p))
            hits +=len(t_index.query(p[start:end]))
        return hits

    @staticmethod
    def approximate_matchs_using_subseq(p, t, n, ival):
        segment_length = int(round(len(p)/(n+1)))
        t_index = SubseqIndex(t,segment_length, ival)
        all_matches = set()
        hits = 0
        for i in range(ival):
            hits += len(t_index.query(p[i:]))
        return hits