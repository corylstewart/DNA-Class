class EditDistance(object):

    def __init__(self, string_one, string_two):
        if len(string_one) < len(string_two):
            string_one, string_two = string_two, string_one
        self.x = string_one
        self.y = string_two
        self.penalty = {'d':1, 'i':1, 's':1}
        self.dp = [[0] * (len(self.x)+1) for i in range(len(self.y)+1)]
        self.edit_distance = None
        
    def make_edit_distance(self):
        for i in range(len(self.dp[0])-1): self.dp[0][i+1] = self.dp[0][i] + 1
        for i in range(len(self.dp)-1): self.dp[i+1][0] = self.dp[i][0] + 1
        for i in range(1, len(self.dp)):
            for j in range(1, len(self.dp[i])):
                self.dp[i][j] = self.calculate_cell(i, j)
        self.edit_distance = self.dp[-1][-1]
        return self.edit_distance

    def calculate_min_within_large_string(self):
        for i in range(len(self.dp)-1): self.dp[i+1][0] = self.dp[i][0] + 1
        for i in range(1, len(self.dp)):
            for j in range(1, len(self.dp[i])):
                self.dp[i][j] = self.calculate_cell(i, j)
        self.edit_distance = min(self.dp[-1])
        return self.edit_distance


    def calculate_cell(self, i, j):
        left = self.dp[i][j-1] + self.penalty['i']
        top = self.dp[i-1][j] + self.penalty['d']
        if self.x[j-1] != self.y[i-1]:
            pen = self.penalty['s']
        else:
            pen = 0
        diag = self.dp[i-1][j-1] + pen
        return min(left, top, diag)


    def print_dp(self):
        for row in self.dp: print row
        print ''


def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]
