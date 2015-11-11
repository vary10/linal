import itertools


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

    def __setitem__(self, key, value):
        self.matrix[key] = value

    def __getitem__(self, key):
        return self.matrix[key]

    def __str__(self):
        res = ""
        for i in range(self.m):
            res += " ".join(map(str, list(map("{:0.4f}".format, self[i])))) + "\n"
        res = res[:-1]
        return res

    def __iter__(self):
        return self.matrix.__iter__()

    def __eq__(self, other):
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
                    print("Умножение не выполнено.\n\
                              Количество столбцов первой матрицы \
                              не равно количеству столбцов второй.")
            else:
                print("Неправвильное")

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
                res[i][j] = (-1) ** (i + j) * tmp.det()
        return res

    def rearrange(self):
        i = self.m - 1
        while i >= 0:
            if self[i].count(0) == self.n:
                k = self.matrix.pop(i)
                self.matrix.append(k)
                i -=1
            else:
                i -= 1

    def gauss(self, m, n):
        try:
            self.rearrange()
            if self[m][n] != 0:
                Matrix.div_line(self[m], self[m][n])
                for i in range(self.m):
                    if i != m:
                        k = self[i][n] / self[m][n]
                        for j in range(n, self.n):
                            self[i][j] -= self[m][j] * k
            return self.gauss(m + 1, n + 1)
        except:
            self.rearrange()
            return self

    def div_line(line, n):
        for i in range(len(line) - 1, -1, -1):
            line[i] /= n
            if line[i] == -0:
                line[i] = 0.0
        return line

    def s_rk(self):
        self.gauss(0, 0)
        for i in range(self.m):
            if self[i].count(0) == self.n:
                return i
        return self.m

    def rk(self):
        return self.copy((0, self.m)).s_rk()

    def decomp(self):
        tmp = self.copy((0, self.m))
        rank = tmp.s_rk()
        C = tmp.copy((0, rank))
        B = self.ccopy((0, rank))
        return (B, C)

    def __pow__(self, other):
        if other == "*":
            return self.transposed()
        elif other == -1:
            return self.inverse()

    def pseudo_inverse(self):
        B, C = self.decomp()
        res = C.transposed() * (C * C.transposed()).inverse()
        res1 = (B.transposed() * B).inverse() * B.transposed()
        return res * res1

# a = Matrix([[1, 1, 1], [2, 2, 2], [3, 3, 3], [1, 2, 15]])
# b = Matrix([[15], [4], [5]])
# print(a.pseudo_inverse())
