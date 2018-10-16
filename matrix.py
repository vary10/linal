import itertools
from copy import deepcopy

class Matrix:

    def __init__(self, matrix=[]):
        self.matrix = matrix
        if matrix:
            self.m = len(matrix)
            self.n = len(matrix[0])
            self.dim = (self.m, self.n)
        else:
            self.m = 0
            self.n = 0
            self.dim = (0, 0)

    def check(self):
        self.matrix = [
            [0 if abs(e) < 1e-7 else float(e) for e in l]
            for l in self.matrix
        ]
    
    def __setitem__(self, key, value):
        self.matrix[key] = value

    def __getitem__(self, key):
        return self.matrix[key]

    def __str__(self):
        self.check()
        res = ""
        for i in range(self.m):
            res += " ".join(
                map(str, list(map("{:0.6f}".format, self[i])))
            ) + "\n"
        res = res[:-1]
        return res

    def __iter__(self):
        return self.matrix.__iter__()

    def __eq__(self, other):
        self.check()
        other.check()
        return (self.matrix == other.matrix)

    def __add__(self, other):
        if type(other) is Matrix:
            if other.dim == self.dim:
                res = Matrix([[None] * self.n for i in range(self.m)])
                for i in range(self.m):
                    for j in range(self.n):
                        res[i][j] = self[i][j] + other[i][j]
                return res

    def __iadd__(self, other):
        if type(other) is Matrix:
            if other.dim == self.dim:
                self = self + other
        return None

    def __mul__(self, other):
        try:
            val = int(other)
            res = Matrix([[None] * self.n for i in range(self.m)])
            for i in range(self.m):
                for j in range(self.n):
                    res[i][j] = other * self[i][j]
            return res
        except:
            if type(other) is Matrix:
                if self.n == other.m:
                    res = Matrix([[0] * other.n for i in range(self.m)])
                    for i in range(self.m):
                        for j in range(other.n):
                            line = self.copy((i, i + 1))
                            column = other.ccopy((j, j + 1)).transposed()
                            for k in range(line.n):
                                res[i][j] += line[0][k] * column[0][k]
                    return res
                else:
                    print("Wrong dims")
            else:
                print("Nor matrix, nor num")

    def __truediv__(self, other):
        return self * (1 / other)

    def copy(self, tup):
        tmp = []
        for p in range(tup[0], tup[1]):
            tmp += [self[p][:]]
        tmp = Matrix(tmp)
        return tmp

    def ccopy(self, tup):
        tmp = []
        for i in range(self.m):
            tmp += [self[i][tup[0]:tup[1]]]
        tmp = Matrix(tmp)
        return tmp

    def transposed(self):
        res = Matrix([[None] * self.m for i in range(self.n)])
        for i in range(self.m):
            for j in range(self.n):
                res[j][i] = self[i][j]
        return res

    def inversions(perm):
        res = 0
        for i in range(len(perm)):
            for j in range(i + 1, len(perm)):
                if perm[i] > perm[j]:
                    res += 1
        return res

    def sq_norm(self):
        return sum([sum(l) for l in (self ** 2).matrix])
    
    def trace(self):
        return sum([self[i][i] for i in range(self.m)])
    
    def det(self):
        if self.m == self.n:
            res = 0
            for perm in itertools.permutations(range(self.n)):
                tmp = (-1) ** Matrix.inversions(perm)
                for i in range(self.n):
                    tmp *= self[i][perm[i]]
                res += tmp
            return res
        else:
            print("Матрица не квадратная.")

    def inverse(self):
        return self.matrix_of_algebraics().transposed() / self.det()

    def matrix_of_algebraics(self):
        res = Matrix([[None] * self.n for i in range(self.m)])
        for i in range(self.m):
            for j in range(self.n):
                tmp = []
                for p in range(self.m):
                    tmp += [self[p][:]]
                    tmp[p].pop(j)
                tmp.pop(i)
                tmp = Matrix(tmp)
                d = tmp.det()
                res[i][j] = (-1) ** (i + j) * d if d != 0 else d
        return res

    def rearrange(self, idxs):
        i = self.m - 1
        while i >= 0:
            if self[i].count(0) == self.n:
                k = self.matrix.pop(i)
                idx = idxs.pop(i)
                self.matrix.append(k)
                idxs.append(idx)
                i -=1
            else:
                i -= 1
        return idxs

    def gauss(self, m=0, n=0, idxs=0):
        tmp = deepcopy(self)
        if idxs == 0:
            idxs = [i for i in range(self.m)]
        try:
            idxs = tmp.rearrange(idxs)
            if tmp[m][n] != 0:
                Matrix.div_line(tmp[m], tmp[m][n])
                for i in range(tmp.m):
                    if i != m:
                        k = tmp[i][n] / tmp[m][n]
                        for j in range(n, tmp.n):
                            tmp[i][j] -= tmp[m][j] * k
            return tmp.gauss(m + 1, n + 1, idxs)
        except:
            idxs = tmp.rearrange(idxs)
            return tmp, idxs

    def div_line(line, n):
        for i in range(len(line) - 1, -1, -1):
            line[i] /= n
            if line[i] == -0:
                line[i] = 0.0
        return line

    def rk(self):
        tmp, idxs = self.gauss()
        for i in range(tmp.m):
            if tmp[i].count(0) == tmp.n:
                return i, idxs
        return tmp.m, idxs

    def decomp(self):
        tmp = deepcopy(self)
        rank, idxs = tmp.rk()
        C = Matrix([self[i] for i in idxs[:rank]])
        B = self * (C.transposed() * (C * C.transposed()).inverse())
        return (B, C)

    def __pow__(self, other):
        if other == 2:
            return Matrix([[e ** 2 for e in l] for l in self.matrix])
        elif other == "*":
            return self.transposed()
        elif other == -1:
            return self.inverse()

    def pseudo_inverse(self):
        B, C = self.decomp()
        res = C.transposed() * (C * C.transposed()).inverse()
        res1 = (B.transposed() * B).inverse() * B.transposed()
        return res * res1