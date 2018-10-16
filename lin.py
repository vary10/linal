from math import cos, acos, pi
from matrix import *

def cross_product(a, b):
    return Matrix([[
        a[0][1] * b[0][2] - a[0][2] * b[0][1],
        a[0][2] * b[0][0] - a[0][0] * b[0][2],
        a[0][0] * b[0][1] - a[0][1] * b[0][0]
    ]])

def eye(n):
    m = [[0] * n for i in range(n)]
    for i in range(n):
        m[i][i] = 1
    return Matrix(m)

def eigen_vals(A):
    q = A.trace() / 3
    t = (A + eye(3) * -q)
    p = (t.sq_norm() / 6) ** 0.5
    B = t / p
    c = B.det() / 2
    phi = acos(c) / 3
    e1 = q + p * 2 * cos(phi)
    e3 = q + p * 2 * cos(phi + 2/3 * pi)
    e2 = A.trace() - e1 - e3
    return e1, e2, e3

def eigen_vec(A, eigen_val):
    B = A + eye(3) * -eigen_val
    m_v = Matrix([[0, 0, 0]])
    for i in range(0, 2):
        for j in range(i + 1, 3):
            v = cross_product(Matrix([B[i]]), Matrix([B[j]]))
            if v.sq_norm() > m_v.sq_norm():
                m_v = Matrix(v.matrix)
    return m_v / m_v.sq_norm() ** 0.5

def SVD(A):
    S = A.transposed() * A
    evals = eigen_vals(S)
    evecs = Matrix([eigen_vec(S, ev).matrix[0] for ev in evals]).transposed()
    sigmas = [e ** 0.5 for e in evals]
    P = A * evecs
    left_vecs = [
        Matrix([P.transposed().matrix[i]]) / sigmas[i]
        for i in range(len(sigmas))
    ]
    left_vecs = Matrix([v.matrix[0] for v in left_vecs]).transposed()
    sigmas = Matrix([[sigmas[0], 0, 0], [0, sigmas[1], 0], [0, 0, sigmas[2]]])
    print("U:\n", left_vecs)
    print("S:\n", sigmas)
    print("V:\n", evecs)
    return left_vecs, sigmas, evecs

def solve(A, B):
    res = A.pseudo_inverse() * B
    for i in range(res.m):
        print(chr(97 + i), "=", "{:0.5f}".format(res[i][0]))
    print()
    return res

def dec(A):
    B, C = A.decomp()
    print("B:", B, "", sep="\n")
    print("C:", C, "", sep="\n")
    print("A+: ", A.pseudo_inverse(), sep="\n")
    print()
    return B, C, A.pseudo_inverse()