class WordSearch(object):

    @staticmethod
    def naive_find_matches(p, t):
        matches = list()
        for i in xrange(len(t)-len(p)+1):
            matched = True
            for j in range(len(p)):
                if p[j] != t[i+j]:
                    matched = False
                    break
            if matched:
                matches.append(i)
        return matches

    @staticmethod
    def naive_find_matches_count(p, t):
        count = 0
        for i in xrange(len(t)-len(p)+1):
            for j in range(len(p)):
                count += 1
                if p[j] != t[i+j]:
                    matched = False
                    break
        return count

    @staticmethod
    def better_find_matches(p, t):
        matches = list()
        nxt = t[0:].find(p)
        while nxt != -1:
            matches.append(nxt)
            i = t[nxt+1:].find(p)
            if i == -1:
                break
            nxt +=  i + 1
        return matches

    @staticmethod
    def find_matches_with_errors(p, t, e):
        matches = list()
        for i in xrange(len(t)-len(p)+1):
            matched = True
            error_count = 0
            for j in range(len(p)):
                if p[j] != t[i+j]:
                    error_count += 1
                    if error_count > e:
                        matched = False
                        break
            if matched:
                matches.append(i)
        return matches

    @staticmethod
    def make_compliment(p):
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        other = ''
        for letter in p:
            other = comp[letter] + other
        return other