from sympy import *
import math

init_printing()


def divisor_generator(n):
    n = abs(n)
    large_divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i * i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor


def divide_poly(a, b):
    if a == 0:
        return 0, 0
    if degree(a, x) == degree(b, x) == 0 or a == 0:
        return Rational(a, b), 0
    else:
        return div(a, b, domain='QQ')


def row_change(a_temp, i, j):
    a_temp[i, :], a_temp[j, :] = a_temp[j, :], a_temp[i, :]


def col_change(a_temp, i, j):
    a_temp[:, i], a_temp[:, j] = a_temp[:, j], a_temp[:, i]


def row_mult(a_temp, i, pol):
    a_temp[i, :] *= pol


def col_mult(a_temp, i, pol):
    a_temp[:, i] *= pol


def row_add(a_temp, i, j, pol=1):
    a_temp[i, :] += a_temp[j, :] * pol


def col_add(a_temp, i, j, pol=1):
    a_temp[:, i] += a_temp[:, j] * pol


def check_submatrix(a_temp, i):
    place = (i, i)
    cur_pol = a_temp[i, i]
    for j in range(i + 1, a_temp.shape[0]):
        for k in range(i + 1, a_temp.shape[0]):
            q, r = divide_poly(a_temp[j, k], a_temp[i, i])
            if r != 0 and degree(r, x) < degree(cur_pol, x):
                cur_pol = r
                place = (j, k)
    if place != (i, i):
        row_add(a_temp, i, place[0])
        return False
    else:
        return True


def row_col_rem(a_temp, q_temp, i):
    for j in range(i + 1, a_temp.shape[0]):
        if a_temp[j, i] != 0:
            q, r = divide_poly(a_temp[j, i], a_temp[i, i])
            if q != 0:
                row_add(a_temp, j, i, -q)
    for j in range(i, a_temp.shape[0]):
        for k in range(i, a_temp.shape[0]):
            a_temp[j, k] = expand(a_temp[j, k])
    for j in range(i + 1, a_temp.shape[0]):
        if a_temp[i, j] != 0:
            q, r = divide_poly(a_temp[i, j], a_temp[i, i])
            if q != 0:
                col_add(a_temp, j, i, -q)
                col_add(q_temp, j, i, -q)
    for j in range(i, a_temp.shape[0]):
        for k in range(i, a_temp.shape[0]):
            a_temp[j, k] = expand(a_temp[j, k])


def row_col_clean(a_temp, q_temp, i):
    cur_pol = a_temp[i, i]
    place = (i, i)
    check = False
    for j in range(i + 1, a_temp.shape[0]):
        if a_temp[j, i] != 0:
            check = True
            if degree(cur_pol, x) > degree(a_temp[j, i], x):
                place = (j, i)
                cur_pol = a_temp[j, i]
    for j in range(i + 1, a_temp.shape[0]):
        if a_temp[i, j] != 0:
            check = True
            if degree(cur_pol, x) > degree(a_temp[i, j], x):
                place = (i, j)
                cur_pol = a_temp[i, j]
    if not check:
        return False
    else:
        if place[0] != i:
            row_change(a_temp, place[0], i)
        else:
            col_change(a_temp, i, place[1])
            col_change(q_temp, i, place[1])
        return True


def canonical_form(a_pol):
    mat_size = a_pol.shape[0]
    a_temp = a_pol
    q_temp = eye(mat_size)
    for i in range(mat_size - 1):
        place = (i, i)
        cur_pol = a_temp[i, i]
        for j in range(i, mat_size):
            for k in range(i, mat_size):
                if degree(a_temp[j, k], x) < degree(cur_pol, x) and a_temp[j, k] != 0 or cur_pol == 0:
                    cur_pol = a_temp[j, k]
                    place = (j, k)
        if cur_pol == 0:
            break
        row_change(a_temp, place[0], i)
        col_change(a_temp, place[1], i)
        col_change(q_temp, place[1], i)
        while True:
            while row_col_clean(a_temp, q_temp, i):
                row_col_rem(a_temp, q_temp, i)
            if check_submatrix(a_temp, i):
                break
    for i in range(mat_size):
        for j in range(mat_size):
            if degree(a_temp[i, j], x) == 0 and a_temp[i, j] != 0:
                a_temp[i, j] = 1
            if degree(a_temp[i, j], x) > 0:
                a_temp[i, j] = factor(a_temp[i, j])
    return a_temp, q_temp


def rational_roots(pol):
    coeffs_list = Poly(pol).all_coeffs()
    count = 0
    while coeffs_list[len(coeffs_list) - 1 - count] == 0:
        count += 1
    last = coeffs_list[-(count + 1)]
    if count != 0:
        roots = {0: count}
    else:
        roots = {}
    for i in divisor_generator(coeffs_list[0]):
        for j in divisor_generator(last):
            cur_pol = pol
            while cur_pol.subs({x: Rational(j, i)}) == 0:
                roots[Rational(j, i)] = roots.get(Rational(j, i), 0) + 1
                cur_pol, r = divide_poly(cur_pol, x - Rational(j, i))
            cur_pol = pol
            while cur_pol.subs({x: Rational(-j, i)}) == 0:
                roots[Rational(-j, i)] = roots.get(Rational(-j, i), 0) + 1
                cur_pol, r = divide_poly(cur_pol, x + Rational(j, i))
    return roots


x = Symbol('x')
# a_mat = Matrix([[1, 0, 1], [0, 1, -1], [-1, -1, 1]])
# a_pol = Matrix([[3, 1, 2], [6, 4, -1], [-1, -1, 1]])
# a_pol = Matrix([[0, 1, -1, 1], [-1, 2, -1, 1], [-1, 1, 1, 0], [-1, 1, 0, 1]])
mat_size = int(input())
a_mat = Matrix([[int(input()) for i in range(mat_size)] for j in range(mat_size)])
mat_size = a_mat.shape[0]
B = eye(mat_size)
a_pol = B * x - a_mat
a_canon, q_first = canonical_form(a_pol)
c, temp_pol = a_canon[mat_size - 1, mat_size - 1].as_content_primitive()
all_roots = rational_roots(temp_pol)
count = 0
for root in all_roots:
    count += all_roots[root]
if count != degree(a_canon[mat_size - 1, mat_size - 1], x):
    print('NO')
    exit()
j_mat = eye(mat_size)
place = 0
for i in all_roots:
    temp = place
    for j in range(all_roots[i]):
        j_mat[place, place] = i
        place += 1
    for j in range(all_roots[i] - 1):
        j_mat[temp, temp + 1] = 1
        temp += 1
for i in range(mat_size - 2, -1, -1):
    if degree(a_canon[i, i], x) > 0:
        c, temp_pol = a_canon[i, i].as_content_primitive()
        all_roots = rational_roots(temp_pol)
        for i in all_roots:
            temp = place
            for j in range(all_roots[i]):
                j_mat[place, place] = i
                place += 1
            for j in range(all_roots[i] - 1):
                j_mat[temp, temp + 1] = 1
j_pol = B * x - j_mat
j_canon, q_second = canonical_form(j_pol)
q_major = q_first * (q_second ** -1)
max_degree = 0
for i in range(mat_size):
    for j in range(mat_size):
        if q_major[i, j] != 0:
            max_degree = max(max_degree, degree(q_major[i, j], x))
t_mat = zeros(mat_size)
for k in range(max_degree + 1):
    cur_mat = zeros(mat_size)
    for i in range(mat_size):
        for j in range(mat_size):
            deg = degree(q_major[i, j], x)
            if deg == 0:
                if k == 0:
                    cur_mat[i, j] = q_major[i, j]
                continue
            if q_major[i, j] != 0 and deg >= k:
                coeffs_list = Poly(q_major[i, j]).all_coeffs()
                cur_mat[i, j] = coeffs_list[deg - k]
    t_mat += cur_mat * j_mat ** k
print(j_mat)
print(t_mat)
print(det(t_mat))
print(a_mat * t_mat - t_mat * j_mat)
# x = Symbol('x')
# mat_size = int(input())
# B = eye(mat_size)
# a_pol = B * x - a_pol
# Q_first = eye(mat_size)
# A_canon, Q_first = canonical_form(a_pol)

# b = a_canon[2,2]
# b = expand(b)
# c, p = b.as_content_primitive()
