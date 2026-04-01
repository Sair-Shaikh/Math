polRing.<x> = PolynomialRing(ZZ)
powRing.<y> = PowerSeriesRing(QQ, 30)

# Load the point counts
ptct_to_key = {}
f = open("Dataset/point_counts_cubic.csv", "rb")
for i in f:
    j = i.split(",")
    k = tuple([int(j[t]) for t in range(1,13)])
    if k in ptct_to_key:
        ptct_to_key[k].append(j[0])
    else:
        ptct_to_key[k] = [j[0]]
f.close()
print "Loaded", len(ptct_to_key.keys()), "distinct point count lists"


# Compute Weil Polynomials
zetas = []
ambiguous = []
for points in ptct_to_key.keys():
    k1 = powRing([0] + [QQ(points[j-1])/j for j in range(1,len(points)+1)])
    halfLqT = (exp(k1)*(1-y)*(1-2*y)*(1-4*y)).inverse()
    m = [2^(1-j)*halfLqT[j] for j in range(len(points)+1)]
    j = 11
    while j < len(m) and m[j]==0: 
        j += 1
    if j >= len(m):
        ambiguous.append([ptct_to_key[points], points])
        continue
    s = sign(m[21-j]) // sign(m[j])
    m1 = m[:]
    for k in range(len(points)+1, 22):
        m1.append(s*m[21-k])
    zetas.append([polRing(m1), [int(j) for j in ptct_to_key[points]]])
print len(zetas), " zeta functions uniquely reconstructed"
print sum(len(elem[0]) for elem in ambiguous), " ambiguous cases remain"



# l4 = []
# l5 = []
# for i in ambiguous:
#     k1 = powRing([0] + [QQ(i[1][j-1])/j for j in range(1,11)])
#     k2 = (exp(k1)*(1-y)*(1-2*y)*(1-4*y)).inverse()
#     m = [2^(1-j)*k2[j] for j in range(11)]
#     m1 = m + m[::-1]
#     m2 = m + [-j for j in m[::-1]]
#     u = [polRing(m1), polRing(m2)]
#     v = []
#     for j in u:
#         ans, count = roots_on_unit_circle(j, n=11)
#         if ans == [j] and j(1).is_square() and ej_test(j):
#             v.append(j)
#     if len(v) == 0: 
#         raise RuntimeError
#     if len(v) == 1:
#         l4.append([v[0], [int(k) for k in i[0]]])
#     if len(v) == 2:
#         print u[0].factor(), u[1].factor(), u[0](1), u[1](1)
#         l5.append(i)
# print len(l4), " zeta functions reconstructed with extra constraints"
# f = open("k3final.txt", "wb")
# for i in zetas+l4:
#     f.write(str(i))
#     f.write("\n")
# f.close()
# print sum(len(i[0]) for i in l5), " more ambiguous cases remain"