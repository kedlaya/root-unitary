load("roots_on_unit_circle.sage")
polRing.<x> = PolynomialRing(Integers())

u = polRing([3, 5, 6, 7, 5, 4, 2, -1, -3, -5, -5, -5, -5, -3, -1, 2, 4, 5, 7, 6, 5, 3])

for i in range(2,6):
    for j in range(1, 6):
        time ans, count = roots_on_unit_circle(u, modulus=3^i, n=j, return_nodes=True)
        print i, j, len(ans), count

u = polRing([2401, -343, -5439, -1050, 7156, 5043, -5829, -7990, 1437, 6348, 2115, -332, -1756, -4639, -1802, 3938, 4762, 16, -3366, -2658, -2051, 1572, 5810, 2097, -5558, -3955, 2598, 1931, -831, 1931, 2598, -3955, -5558, 2097, 5810, 1572, -2051, -2658, -3366, 16, 4762, 3938, -1802, -4639, -1756, -332, 2115, 6348, 1437, -7990, -5829, 5043, 7156, -1050, -5439, -343, 2401])

time ans, count = roots_on_unit_circle(u, modulus=7^2, n=28, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^3, n=25, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^3, n=24, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^4, n=20, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^4, n=19, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^4, n=18, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^4, n=17, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^4, n=16, answer_count=2, return_nodes=True)
print len(ans), count
time ans, count = roots_on_unit_circle(u, modulus=7^5, n=1, answer_count=2, return_nodes=True)
print len(ans), count
