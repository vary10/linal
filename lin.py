from matrix import Matrix
from math import *

def solve(A, B, n):
    res = A.pseudo_inverse() * B
    print("задача", n)
    for i in range(res.m):
        print(chr(97 + i), "=", "{:0.4f}".format(res[i][0]))
    print()

def dec(A, n):
    B, C = A.decomp()
    print("задача", n)
    print("B:", B, "", sep="\n")
    print("C:", C, "", sep="\n")
    print("A+: ", A.pseudo_inverse(), sep="\n")
    print()

A2 = Matrix([[1, 1], [3, 1], [5, 1]])
B2 = Matrix([[15], [4], [5]])
solve(A2, B2, 2)


A3 = Matrix([[1, sin(1)], [3, sin(3)], [5, sin(5)], [15, sin(15)]])
B3 = Matrix([[2], [4], [5], [0]])
solve(A3, B3, 3)

A8 = Matrix([[2, 3, 2,], [3, 4, -1]])
B8 = Matrix([[7], [6]])
solve(A8, B8, 8)

A9 = Matrix([[1, -3, 0, 1], [0, 2, -3, 0], [1, -2, 1, 1], [1, 0, -2, 1]])
B9 = Matrix([[-1], [-1], [0], [15]])
solve(A9, B9, 9)

A5a = Matrix([[1, 0, 15]])
dec(A5a, "5a")

A5b = Matrix([[0], [1], [15], [3]])
dec(A5b, "5b")

A5c = Matrix([[1, 0], [-1, 0], [-15, 0], [2, 1]])
dec(A5c, "5c")

A6 = Matrix([[1, 1, 1], [2, 2, 2], [3, 3, 3], [1, 2, 15]])
dec(A6, "6")

f = open("in.txt", "r")
temp = Matrix([list(map(int, f.read().split()))])
temp = temp.transposed()

years = list(range(1949, 2013))
years = [1937] + years
years.pop(years.index(1956))
years.pop(years.index(1971))
years.pop(years.index(1972))
years.pop(years.index(1987))
years.pop(years.index(1958))

print("задача 4:")

const = Matrix([[1] for i in range(60)])
print("если климат стабилен:", const.pseudo_inverse() * temp)

lin = Matrix([years]).transposed().matrix
for i in range(len(lin)):
    lin[i] += [1]
lin = Matrix(lin)
res = lin.pseudo_inverse() * temp
print("если учитывать потепление:", "{:0.4f}".format(res[0][0] * 2015 + res[1][0]))
