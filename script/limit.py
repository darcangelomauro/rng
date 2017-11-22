import math
from numpy import matrix
from numpy import linalg
import numpy
import copy
import itertools


# powers gives ls^n, where ls is a list of characters
def powers(ls, n):
    ls1 = []
    for i in range(n):
        ls1.append(ls)

    res = itertools.product(*ls1)
    res1 = []
    for comb in res:
        s = ''
        for item in comb:
            s += item
        res1.append(s)
    return res1

# derivative gives a list of the form [[x1,y1], [x2, y2], ... ]
# where xi is the polynomial and yi is the coefficient
def derivative(idx, ls):
    ls_der = []

    for item in ls:
        full = []
        j = []
        n = len(item)
        for i in range(n):
            full.append(item[i])
            if item[i] == idx:
                j.append(i)
        # now I know that the list full has 'idx' in position j[0], j[1], ...
        for index in j:
            temp = full[:]
            temp.pop(index)
            ls_temp = []
            ls_temp.append(''.join(temp))
            ls_temp.append(item)
            ls_der.append(ls_temp)

    return ls_der

def slice2D(idx1, idx2, ls_der):

    ls_out = []

    for item in ls_der:
        check = 1
        for i in range(len(item[0])):
            if item[0][i] != idx1 and item[0][i] != idx2:
                check = 0
        if check:
            ls_out.append(item)

    return ls_out

# returns the first (alphabetically ordered) cyclic permutation of elements of a
def cyclic_sort(a):
    
    n = len(a)
    b = [[a[i - j] for i in range(n)] for j in range(n)]
    b = sorted(b)

    return b[0]

# rewrites the elements of a by imposing commutativity of a[i][0] and
# cyclic symmetry of a[i][1]
def rewrite(a):

    b = []

    for item in a:
        com = [item[0][i] for i in range(len(item[0]))]
        cyc = [item[1][i] for i in range(len(item[1]))]

        com = sorted(com)
        cyc = cyclic_sort(cyc)
        
        temp = ["".join(com), "".join(cyc)]
        b.append(temp)

    return b

def second(a):
    return a[1]

def gamma_check(gamma, a):
    
    b = []

    for item in a:
        spl = [item[1][i] for i in range(len(item[1]))]
        mat = numpy.identity( (gamma[0].shape)[0] )

        for idx in spl:
            mat = numpy.dot(mat, gamma[int(idx)])

        tr = numpy.trace(mat)

        if tr != 0:
            b.append([item[0], item[1], str(tr)])

    return b




def shrink(a):
    
    b = []
    for src in a:

        # check if element has already been accounted for
        check = 1
        for item in b:
            if src[1] == item[1]:
                check = 0
        # if not, find the sum
        if check:
            count = 0
            for comp in a:
                if src[1] == comp[1]:
                    count += 1
            b.append([str(count) + " * " + src[0], src[1]])

    return b


# given an alphabet and two indices, gives the 2d slices corresponding to those coordinates 
def masterchef(alphabet, gamma, idx1, idx2):

    ls2 = powers(ls, 2)
    ls4 = powers(ls, 4)
    
    ls2_der = [derivative(str(i), ls2) for i in range(len(ls))]
    ls4_der = [derivative(str(i), ls4) for i in range(len(ls))]

    #ls2_slice = [slice2D(str(idx1), str(idx2), ls2_der[i]) for i in range(len(ls))]
    #ls4_slice = [slice2D(str(idx1), str(idx2), ls4_der[i]) for i in range(len(ls))]
    ls2_slice = [ls2_der[i] for i in range(len(ls))]
    ls4_slice = [ls4_der[i] for i in range(len(ls))]

    ls2_slice = [rewrite(item) for item in ls2_slice]
    ls4_slice = [rewrite(item) for item in ls4_slice]

    ls2_slice = [sorted(item, key = second) for item in ls2_slice]
    ls4_slice = [sorted(item, key = second) for item in ls4_slice]

    ls2_slice = [shrink(item) for item in ls2_slice]
    ls4_slice = [shrink(item) for item in ls4_slice]
    
    ls2_slice = [gamma_check(gamma, item) for item in ls2_slice]
    ls4_slice = [gamma_check(gamma, item) for item in ls4_slice]

    return [ls2_slice, ls4_slice]

# (3,1) geometry
'''
g0 = matrix([[1.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0],[0.0, 0.0, -1.0, -0.0],[0.0, 0.0, -0.0, -1.0]])
g1 = matrix([[0.0, 0.0, 1.0, 0.0],[0.0, 0.0, 0.0, 1.0],[1.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0]])
g2 = matrix([[0j, 0j, 1j, 0j],[0j, (-0+0j), 0j, (-0-1j)],[-1j, 0j, 0j, 0j],[0j, 1j, 0j, (-0+0j)]])
g3 = matrix([[0j, -1j, 0j, 0j],[1j, 0j, 0j, 0j],[0j, 0j, 0j, -1j],[0j, 0j, 1j, 0j]])
g4 = matrix([[0j, (1+0j), 0j, 0j],[(1+0j), 0j, 0j, 0j],[0j, 0j, 0j, (-1+0j)],[0j, 0j, (-1+0j), 0j]])
g5 = matrix([[0j, 0j, 0j, (1+0j)],[0j, 0j, (1+0j), 0j],[0j, (1+0j), 0j, 0j],[(1+0j), 0j, 0j, 0j]])

gamma = [g0, g1, g2, g3, g4, g5]

ls = ['0', '1', '2', '3', '4', '5']
'''


# (2,2) geometry
#'''
g0 = matrix([[1.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0],[0.0, 0.0, -1.0, -0.0],[0.0, 0.0, -0.0, -1.0]])
g1 = matrix([[0.0, 0.0, 1.0, 0.0],[0.0, 0.0, 0.0, 1.0],[1.0, 0.0, 0.0, 0.0],[0.0, 1.0, 0.0, 0.0]])
g2 = matrix([[(1+0j), 0j, 0j, 0j],[0j, (-1+0j), 0j, 0j],[0j, 0j, (1+0j), 0j],[0j, 0j, 0j, (-1+0j)]])
g3 = matrix([[0j, -1j, 0j, 0j],[1j, 0j, 0j, 0j],[0j, 0j, 0j, -1j],[0j, 0j, 1j, 0j]])
g4 = matrix([[0j, 0j, (-1+0j), 0j],[0j, -0j, 0j, (1-0j)],[(1+0j), 0j, 0j, 0j],[0j, (-1+0j), 0j, -0j]])
g5 = matrix([[0j, 0j, 0j, 1j],[(-0+0j), 0j, (-0-1j), 0j],[0j, -1j, 0j, 0j],[1j, 0j, (-0+0j), 0j]])

gamma = [g0, g1, g2, g3]

ls = ['0', '1', '2', '3']
#'''

equations = masterchef(ls, gamma, 2, 3)

for power in equations:
    print("degree: " + str(2*(equations.index(power) + 1)) )
    for item in power:
        if item != []:
            print("eq. number: " + str(power.index(item)) )
            print(item)











