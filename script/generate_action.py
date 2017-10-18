import math
from numpy import matrix
from numpy import linalg
import numpy
import copy
import itertools


class block_trace:
    "block_trace is a class which stores the three trace terms of a block"

    def __init__(self, tr_gamma, tr_im, tr_arg1, tr_arg2, name):
        self.tr_gamma = tr_gamma
        self.tr_im = tr_im
        self.tr_arg1 = tr_arg1
        self.tr_arg2 = tr_arg2
        self.name = name
        if type(tr_gamma) != int:
            print("block_trace error: tr_gamma is not an int")
        if type(tr_im) != str:
            print("block_trace error: tr_im is not a string")
        if type(tr_arg1) != str:
            print("block_trace error: tr_arg1 is not a string")
        if type(tr_arg2) != str:
            print("block_trace error: tr_arg2 is not a string")
        if type(name) != str:
            print("block error: name is not a string")

    def trace_string_format(self):
        return str(self.tr_gamma) + self.tr_im + "*" + self.tr_arg1 + "*" + self.tr_arg2


    def literal_trace_string_format(self):
        if self.tr_im != "":
            return "IU" + "*" + self.tr_arg1 + "*" + self.tr_arg2
        else:
            return self.tr_arg1 + "*" + self.tr_arg2



class block:
    "a block is an object of the type (GAMMA x H x Id) where x denotes tensor product"

    def __init__(self, gamma, arg1, arg2, sign, im, name):
        self.gamma = gamma # this is the gamma matrix
        self.arg1 = arg1 # this is H or L matrix or identity
        self.arg2 = arg2 # same as arg1
        self.sign = sign # this is an overall sign
        self.im = im # if this is 1, there's an overall i factor
        self.name = name # this is the name of the block, just for convenience
        if type(gamma) != matrix:
            print("block error: gamma is not a matrix")
        if type(arg1) != str:
            print("block error: arg1 is not a string")
        if type(arg2) != str:
            print("block error: arg2 is not a string")
        if type(name) != str:
            print("block error: name is not a string")
        if type(sign) != int:
            print("block error: sign is not an int")
        if abs(sign) != 1:
            print("block error: sign doesn't have magnitude 1")
        if type(im) != int:
            print("block error: im is not an int")



    def create_block_trace(self):
        tr = []
        f = ""
        if self.im:
            a = self.sign*self.gamma.trace()*1j
        else:
            a = self.sign*self.gamma.trace()
        if numpy.real(a) != 0 and numpy.imag(a) != 0:
            print("error, trace of gamma matrix is not purely real or purely imaginary")
        b = ""
        c = ""



        if(numpy.imag(a) != 0):
            f = "i"
            a = int(numpy.imag(a))
        else:
            f = ""
            a = int(numpy.real(a))

        if check_arg_zero(self.arg1, "_"):
            a = 0
        elif self.arg1 == "Id":
            b = "dim"
        else:
            b = "tr" + order_arg(self.arg1, "_")

        if check_arg_zero(self.arg2, "_"):
            a = 0
        elif self.arg2 == "Id":
            c = "dim"
        else:
            c = "tr" + order_arg(self.arg2, "_")

        if a == 0:
            tr.append(a)
            tr.append(str(a))
            tr.append(str(a))
            tr.append(str(a))
        else:
            tr.append(a)
            tr.append(f)
            tr.append(b)
            tr.append(c)

        tr = tr[0:2] + sorted(tr[2:])
        return block_trace(tr[0], tr[1], tr[2], tr[3], self.name)


    def block_print(self):
        print(self.sign, end=" ")
        if(self.im):
            print("i", end=" ")
        else:
            print("", end= " ")
        print(self.gamma, end=" ")
        print(self.arg1, end=" ")
        print(self.arg2)



def trace_cleanup(a):
    if type(a) != list:
        print("error, trying to clean up something which is not a list")
    for item in a:
        if type(item) != block_trace:
            print("error, items to clean up are not block_trace")

    literal = []
    gamma = []
    cumul = []
    for item in a:
        literal.append(item.literal_trace_string_format())
        gamma.append(item.tr_gamma)
        cumul.append(0)

    check = 0
    for i in range(0, len(literal)):
        for j in range(0, len(literal)):
            if literal[i] == literal[j]:
                check = 1
                cumul[i] += gamma[j]
        if check == 0:
            print("error, somethings off")

    for i in range(0, len(literal)):
        if cumul[i] == 0:
            literal[i] = "0*0"

    clean = []
    for i in range(0, len(literal)):
        clean.append(str(cumul[i]) + "*" + literal[i])

    clean = [item for item in clean if item != "0*0*0"]

    return list(set(clean))






# computes product of two blocks
def block_product(a, b):
    if type(a) != block or type(b) != block:
        print("error, attempted block product between two non-block objects")
        return
    else:
        g = numpy.matmul(a.gamma, b.gamma)
        s = a.sign*b.sign
        n = a.name + b.name
        arg1 = ""
        arg2 = ""
        if a.arg1 == "Id":
            arg1 = b.arg1
        elif b.arg1 == "Id":
            arg1 = a.arg1
        else:
            arg1 = a.arg1 + "_" + b.arg1
        if a.arg2 == "Id":
            arg2 = b.arg2
        elif b.arg2 == "Id":
            arg2 = a.arg2
        else:
            arg2 = a.arg2 + "_" + b.arg2
        im = 0
        if a.im and b.im:
            s *= -1
        elif a.im or b.im:
            im = 1

        return block(g, arg1, arg2, s, im, n)

def block_multi_product(block_list):
    if type(block_list) != list:
        print("error, attempted multi product with non-list object")
        return
    for item in block_list:
        if type(item) != block:
            print("error, attempted multi product with list of non-block objects")
            return

    n = len(block_list)
    if n < 2:
        print("error, give me at least 2 objects in the multi product")
        return

    res = block_product(block_list[0], block_list[1])
    for i in range(2, n):
        res = block_product(res, block_list[i])

    return res




# takes as input a list of blocks [B1, B2, ..., Bk]
# which represents the dirac operator: B1 + B2 + ... + Bk
# and computes the n-th power of it
def dirac_power(a, n):
    if type(a) != list:
        print("error, attempted power with non-list object")
        return
    for item in a:
        if type(item) != block:
            print("error, attempted power with list of non-block objects")
            return
    if n < 2:
        print("error, attempted power with exponent < 2")
        return

    # this is to generate every combination of addends in n-th power
    res2 = generate_combinations(a, n)
    # now res2 is a list of blocks that have to be multiplied

    # the following performs the block multiplication
    res = []
    for term in res2:
        res.append(block_multi_product(term))

    return res

def generate_combinations(a, n):
    A = []
    for i in range(n):
        A.append(a)
    res1 = itertools.product(*A)
    res2 = [list(item) for item in res1]
    return res2


def order_arg(a, separator):
    if type(separator) != str:
        print("error, separator is not a string")
        return

    s_a = a.split(separator)
    s_a = [item for item in s_a if item != ""]

    # now account for reverse order if a transpose is found
    if a.find("T") != -1:
        s_a_rev = [item for item in reversed(s_a)]
        s_a = [item.replace("T", "") for item in s_a_rev]

    n = len(s_a)
    b = [[s_a[i - j] for i in range(n)] for j in range(n)]
    b.sort()

    c = ""
    for i in range(n):
        c += "." + b[0][i]

    return c

def check_arg_zero(a, separator):
    if type(separator) != str:
        print("error, separator is not a string")
        return

    s_a = a.split(separator)
    s_a = [item for item in s_a if item != ""]

    n = len(s_a)
    if n == 1:
        if s_a[0].find("L") != -1:
            return 1

    return 0

# returns the structure of matrix products
# takes as input the action (as given by trace_cleanup) and list of existing products
def product_structure(action, existing):
    # isolate matrix products in the action
    s_action = [item.split("*") for item in action]
    s_action = [[item for item in item1 if item.find("tr") != -1] for item1 in s_action]
    s_s_action = [[item.split(".") for item in item1] for item1 in s_action]
    s_s_action = [[["." + item for item in item1 if item != "tr"] for item1 in item2] for item2 in s_s_action]

    # now s_s_action is a list of the form [ [[.Hi, .Lj]], [[.Hk], [.Lm]], ... ]
    # the biggest list is the whole action, the middle list is an addend, the
    # smallest list is a trace (so in the example it would be: trHiLj + trHk*trLm)

    # non_existing will be the list of yet non-existing matrix products
    non_existing = []
    for addend in s_s_action:
        for trace in addend:
            name = ""
            for matrix in trace:
                name += matrix
            non_existing.append(name)

    # remove duplicates and order non_existing based on length (from shorter to longer)
    non_existing = list(set(non_existing))
    non_existing.sort(key = lambda s: len(s))

    # now non_existing is an ordered list of the form [.Hn, .Ly, .Hi.Lj, .Hk.Hm, ..... ]
    # which holds every matrix product we need in order of length

    # now split products into two products in such a way as to obtain all the necessary products
    # and repeat until every product has been split in a product of two elementary matrices
    stop_check = 1
    while stop_check:
        non_existing_copy = [product.split(".") for product in non_existing]
        for product in non_existing_copy:
            product.remove("")
        non_existing_s = []
        for product in non_existing_copy:
            n = int(len(product)/2.)
            pair1 = []
            pair2 = []
            for i in range(0,n):
                pair1.append(product[i])
            for i in range(n, len(product)):
                pair2.append(product[i])
            non_existing_s.append([["." + ".".join(pair1)], ["." + ".".join(pair2)]])
        non_existing_s = [[pair for pair in product if pair != ["."]] for product in non_existing_s]

        # now non_existing_s is a list of the form [[[.Hi], [.LjLk]], ... ] which holds every necessary
        # matrix product in the action. but the list structure is kinda redundant, so now we flatten

        necessary = []
        for product in non_existing_s:
            for pair in product:
                for matrix in pair:
                    necessary.append(matrix)
        necessary = list(set(necessary))
        necessary.sort(key = lambda s: len(s))

        """
        for item in necessary:
            if not item in non_existing:
                non_existing.append(item)
        non_existing.sort(key = lambda s: len(s))
        """
        stop_check = 0
        for item in necessary:
            if item.count(".") > 2:
                stop_check = 1

    necessary = [item for item in necessary if not item in non_existing]
    # so now the list necessary holds every ghost matrix product we need to compute
    # the action (ghost means we don't need their trace)

    # structure will be the list holding information on how to build the matrices in existing
    structure = []
    for item in non_existing:
        for item1 in existing:
            if item == item1:
                structure.append([item1])

    # lcurrent is the maximum number of matrices multiplied together currently existing
    lcurrent = 0
    for item in existing:
        count = item.count(".")
        if count > lcurrent:
            lcurrent = count

    # now we count the maximum number of matrices to be multiplied together: lmax
    lmax = 0
    for product in non_existing:
        count = product.count(".")
        if count > lmax:
            lmax = count

    # stop the cycle when lcurrent == lmax
    ghost = []
    while lcurrent < lmax:
        combinations = generate_combinations(existing, 2)
        for item in non_existing:
            for item1 in combinations:
                if item == "".join(item1):
                    structure.append(item1)
                    existing.append("".join(item1))
                    ghost.append(0)
                    break
        for item in necessary:
            for item1 in combinations:
                if item == "".join(item1):
                    structure.append(item1)
                    existing.append("".join(item1))
                    ghost.append(1)
                    break
        non_existing = [item for item in non_existing if not item in existing]
        necessary = [item for item in necessary if not item in existing]
        for item in existing:
            count = item.count(".")
            if count > lcurrent:
                lcurrent = count

    if len(non_existing) != 0:
        print("error, still non-empty non-existing list:")
        print(non_existing)
        return


    structure = [item for item in structure if len(item)>1]

    existing = [item for item in existing if item.count(".") > 1]
    existing1 = []
    structure1 = []
    hermitian = []
    for item in existing:
        existing1.append(item.replace(".", ""))
    for item in structure:
        structure1.append([item1.replace(".", "") for item1 in item])
    for item in structure:
        if check_gh(item[0]):
            hermitian.append("hleft")
        elif check_gh(item[1]):
            hermitian.append("hright")
        else:
            hermitian.append("g")
    #print(existing1)
    #print(structure1)

    for i in range(len(structure1)):
        structure1[i].insert(0, existing1[i])
        structure1[i].insert(len(structure1[i]), hermitian[i])
        structure1[i].insert(len(structure1[i]), ghost[i])
    return structure1




def needed(comb, non_existing, existing):
    for item in existing:
        if comb == item:
            return 0

    for item in non_existing:
        temp = item.split(".")
        temp = [item1 for item1 in temp if item1 != ""]
        pairs = []
        for i in range(int(len(temp)/2.)):
            pairs.append(item[2*i] + "." + item[2*i+1] + ".")
        for item1 in pairs:
            if comb == item1:
                return 1

    return 0





# check_gh checks whether to use hermitian multiplication or not
def check_gh(prod):
    if type(prod) != str:
        print("error, check_gh need a string object to evaluate")
        return

    terms = prod.split(".")
    terms = [item for item in terms if item != ""]
    rev_terms = [item for item in reversed(terms)]

    check = 1

    for i in range(len(terms)):
        if terms[i] != rev_terms[i]:
            check = 0

    return check


def REAL_action(cl):
    if type(cl) != list:
        print("error, the action must be in list form")

    action = ""

    for item in cl:
        if item.count("tr") == 1:
            item_s = item.split("tr")
            imag_check = 0
            if item_s[0].find("IU") != -1:
                imag_check = 1
                item_s[0] = item_s[0].replace("IU*", "")
            item_s[1] = item_s[1].replace(".", "")
            if item_s[0][0] == "-":
                if imag_check:
                    action += item_s[0] + "GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), tr" + item_s[1] + "))"
                else:
                    action += item_s[0] + "GSL_REAL(tr" + item_s[1] + ")"
            else:
                if imag_check:
                    action += "+" + item_s[0] + "GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.), tr" + item_s[1] + "))"
                else:
                    action += "+" + item_s[0] + "GSL_REAL(tr" + item_s[1] + ")"
        elif item.count("tr") == 2:
            item_s = item.split("tr")
            imag_check = 0
            if item_s[0].find("IU") != -1:
                imag_check = 1
                item_s[0] = item_s[0].replace("IU*", "")
            item_s[1] = item_s[1].replace("*", "")
            item_s[1] = item_s[1].replace(".", "")
            item_s[2] = item_s[2].replace(".", "")
            if item_s[0][0] == "-":
                if imag_check:
                    action += item_s[0] + "GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + ")))"
                else:
                    action += item_s[0] + "GSL_REAL(gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + "))"
            else:
                if imag_check:
                    action += "+" + item_s[0] + "GSL_REAL(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + ")))"
                else:
                    action += "+" + item_s[0] + "GSL_REAL(gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + "))"
    return action

def IMAG_action(cl):
    if type(cl) != list:
        print("error, the action must be in list form")

    action = ""

    for item in cl:
        if item.count("tr") == 1:
            item_s = item.split("tr")
            imag_check = 0
            if item_s[0].find("IU") != -1:
                imag_check = 1
                item_s[0] = item_s[0].replace("IU*", "")
            item_s[1] = item_s[1].replace(".", "")
            if item_s[0][0] == "-":
                if imag_check:
                    action += item_s[0] + "GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), tr" + item_s[1] + "))"
                else:
                    action += item_s[0] + "GSL_IMAG(tr" + item_s[1] + ")"
            else:
                if imag_check:
                    action += "+" + item_s[0] + "GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.), tr" + item_s[1] + "))"
                else:
                    action += "+" + item_s[0] + "GSL_IMAG(tr" + item_s[1] + ")"
        elif item.count("tr") == 2:
            item_s = item.split("tr")
            imag_check = 0
            if item_s[0].find("IU") != -1:
                imag_check = 1
                item_s[0] = item_s[0].replace("IU*", "")
            item_s[1] = item_s[1].replace("*", "")
            item_s[1] = item_s[1].replace(".", "")
            item_s[2] = item_s[2].replace(".", "")
            if item_s[0][0] == "-":
                if imag_check:
                    action += item_s[0] + "GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + ")))"
                else:
                    action += item_s[0] + "GSL_IMAG(gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + "))"
            else:
                if imag_check:
                    action += "+" + item_s[0] + "GSL_IMAG(gsl_complex_mul(gsl_complex_rect(0., 1.),gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + ")))"
                else:
                    action += "+" + item_s[0] + "GSL_IMAG(gsl_complex_mul(tr" + item_s[1] + ",tr" + item_s[2] + "))"
    return action




# DPL is a list that must be the output of the dirac_power function
# init is a list of the initial global H and L variables
# (so they are an array of matrices), and init_mod is a list of the
# initial H and L variables but without square brackets (that can cause problem otherwise)
# if print_check, then it prints out the whole variable declarations, otherwise it just returns the
# real and imaginary part (as a list, [RE, IM]) of the trace as a string in a c-like return statement
def masterchef(DPL, init, init_mod, print_check):
    ls = []
    for item in DPL:
        ls.append(item.create_block_trace())

    cl = trace_cleanup(ls)
    print(cl)

    init_mod2 = ["." + item for item in init_mod]
    structure = product_structure(cl, init_mod2)
    #print(structure)

    if print_check:
        print("")
        print("// HERE IT STARTS ****************")
        print("")
        for item in init_mod:
            print("gsl_matrix_complex* " + item + " = " + init[init_mod.index(item)] + ";")
            if item.find("L") == -1:
                print("gsl_complex tr" + item + " = trace(" + item + ");")
        for item in structure:
            print("gsl_matrix_complex* " + item[0] + " = gsl_matrix_complex_calloc(dim, dim);")
            if item[3] == "hleft":
                print("gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, " + item[1] + ", " + item[2] + ", GSL_COMPLEX_ZERO, " + item[0] + ");")
            elif item[3] == "hright":
                print("gsl_blas_zhemm(CblasRight, CblasUpper, GSL_COMPLEX_ONE, " + item[2] + ", " + item[1] + ", GSL_COMPLEX_ZERO, " + item[0] + ");")
            else:
                print("gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, " + item[1] + ", " + item[2] + ", GSL_COMPLEX_ZERO, " + item[0] + ");")
            if not item[4]:
                print("gsl_complex tr" + item[0] + " = trace(" + item[0] + ");")
        for item in structure:
            print("gsl_matrix_complex_free(" + item[0] + ");")

    return [REAL_action(cl), IMAG_action(cl)]


"""
# (0,3) geom
init = ["H[0]", "L[0]", "L[1]", "L[2]"]
init_mod = ["H0", "L1", "L2", "L3"]
g0 = matrix([[1,0], [0,1]])
g1 = matrix([[1j,0], [0,-1j]])
g2 = matrix([[0,1], [-1,0]])
g3 = matrix([[0,1j], [1j,0]])
a0 = block(g0, "H0", "Id", 1, 0, "a0")
a1 = block(g0, "Id", "TH0", 1, 0, "a1")
b0 = block(g1, "L1", "Id", 1, 1, "b0")
b1 = block(g1, "Id", "TL1", -1, 1, "b1")
c0 = block(g2, "L2", "Id", 1, 1, "c0")
c1 = block(g2, "Id", "TL2", -1, 1, "c1")
d0 = block(g3, "L3", "Id", 1, 1, "d0")
d1 = block(g3, "Id", "TL3", -1, 1, "d1")

list1 = [a0, a1, b0, b1, c0, c1, d0, d1]
"""


"""
# (0,2) geom
init = ["L[0]", "L[1]"]
init_mod = ["L0", "L1"]
g0 = matrix([[1j,0], [0,-1j]])
g1 = matrix([[0,1], [-1,0]])
a0 = block(g0, "L0", "Id", 1, 1, "a0")
a1 = block(g0, "Id", "L0", -1, 1, "a1")
b0 = block(g1, "L1", "Id", 1, 1, "b0")
b1 = block(g1, "Id", "L1", -1, 1, "b1")

list1 = [a0, a1, b0, b1]
"""

#"""
# (2,0) geom
init = ["H[0]", "H[1]"]
init_mod = ["H0", "H1"]
g0 = matrix([[1,0], [0,-1]])
g1 = matrix([[0,1], [1,0]])
a0 = block(g0, "H0", "Id", 1, 0, "a0")
a1 = block(g0, "Id", "TH0", 1, 0, "a1")
b0 = block(g1, "H1", "Id", 1, 0, "b0")
b1 = block(g1, "Id", "TH1", 1, 0, "b1")

list1 = [a0, a1, b0, b1]
#"""

"""
# (1,1) geom
init = ["H[0]", "L[0]"]
init_mod = ["H0", "L0"]
g0 = matrix([[1,0], [0,-1]])
g1 = matrix([[0,1], [-1,0]])
a0 = block(g0, "H0", "Id", 1, 0, "a0")
a1 = block(g0, "Id", "TH0", 1, 0, "a1")
b0 = block(g1, "L0", "Id", 1, 1, "b0")
b1 = block(g1, "Id", "TL0", -1, 1, "b1")

list1 = [a0, a1, b0, b1]
"""

"""
# (0,2) geom
init = ["L[0]", "L[1]"]
init_mod = ["L0", "L1"]
g0 = matrix([[1j,0], [0,-1j]])
g1 = matrix([[0,1], [-1,0]])
a0 = block(g0, "L0", "Id", 1, 1, "a0")
a1 = block(g0, "Id", "TL0", -1, 1, "a1")
b0 = block(g1, "L1", "Id", 1, 1, "b0")
b1 = block(g1, "Id", "TL1", -1, 1, "b1")

list1 = [a0, a1, b0, b1]
"""



"""
# (1,3) geom
init = ["H[0]", "L[0]", "L[1]", "L[2]", "L[3]", "L[4]", "L[5]", "H[1]"]
init_mod = ["H0", "L0", "L1", "L2", "L3", "L4", "L5", "H1"]

g0 = matrix([[0,0,1,0], [0,0,0,1], [1,0,0,0], [0,1,0,0]])
g1 = matrix([[0,0,0,-1], [0,0,-1,0], [0,1,0,0], [1,0,0,0]])
g2 = matrix([[0,0,0,1j], [0,0,-1j,0], [0,-1j,0,0], [1j,0,0,0]])
g3 = matrix([[0,0,-1,0], [0,0,0,1], [1,0,0,0], [0,-1,0,0]])
g4 = matrix([[0,0,-1j,0], [0,0,0,1j], [-1j,0,0,0], [0,1j,0,0]])
g5 = matrix([[0,0,0,1], [0,0,-1,0], [0,1,0,0], [-1,0,0,0]])
g6 = matrix([[0,0,0,-1j], [0,0,-1j,0], [0,-1j,0,0], [-1j,0,0,0]])
g7 = matrix([[0,0,1j,0], [0,0,0,1j], [-1j,0,0,0], [0,-1j,0,0]])

a0 = block(g0, "H0", "Id", 1, 0, "a0")
a1 = block(g0, "Id", "TH0", 1, 0, "a1")
b0 = block(g1, "L0", "Id", 1, 1, "b0")
b1 = block(g1, "Id", "TL0", -1, 1, "b1")
c0 = block(g2, "L1", "Id", 1, 1, "c0")
c1 = block(g2, "Id", "TL1", -1, 1, "c1")
d0 = block(g3, "L2", "Id", 1, 1, "d0")
d1 = block(g3, "Id", "TL2", -1, 1, "d1")
e0 = block(g4, "L3", "Id", 1, 1, "e0")
e1 = block(g4, "Id", "TL3", -1, 1, "e1")
f0 = block(g5, "L4", "Id", 1, 1, "f0")
f1 = block(g5, "Id", "TL4", -1, 1, "f1")
h0 = block(g6, "L5", "Id", 1, 1, "h0")
h1 = block(g6, "Id", "TL5", -1, 1, "h1")
i0 = block(g7, "H1", "Id", 1, 0, "i0")
i1 = block(g7, "Id", "TH1", 1, 0, "i1")
list1 = [a0, a1, b0, b1, c0, c1, d0, d1, e0, e1, f0, f1, h0, h1, i0, i1]
"""

"""
# (0,1) geom
init = ["L[0]"]
init_mod = ["L0"]
g0 = matrix([[-1j]])
a0 = block(g0, "L0", "Id", 1, 1, "a0")
a1 = block(g0, "Id", "TL0", -1, 1, "a1")

list1 = [a0, a1]
"""

"""
# (1,0) geom
init = ["H[0]"]
init_mod = ["H0"]
g0 = matrix([[1]])
a0 = block(g0, "H0", "Id", 1, 0, "a0")
a1 = block(g0, "Id", "TH0", 1, 0, "a1")

list1 = [a0, a1]
"""

"""
#for debug
#for D4+gD2

D4 = dirac_power(list1, 4)
D2 = dirac_power(list1, 2)

trD2 = masterchef(D2, init, init_mod, 0)
trD4 = masterchef(D4, init, init_mod, 1)


print("if(control[0])")
print("vecS[0] = G*(" + trD2[0] + ") " + trD4[0] + ";")
print("if(control[1])")
print("vecS[1] = G*(" + trD2[1] + ") " + trD4[1] + ";")
print("if(control[2])")
print("vecS[2] = GSL_REAL(trH0);")
print("if(control[3])")
print("vecS[3] = GSL_REAL(trH0H0);")
print("if(control[4])")
print("vecS[4] = GSL_REAL(trH1);")
print("if(control[5])")
print("vecS[5] = GSL_REAL(trH1H1);")
print("}")



#for D2
D2 = dirac_power(list1, 2)

trD2 = masterchef(D2, init, init_mod, 1)

print("if(control[0])")
print("vecS[0] = " + trD2[0] + ";")
print("if(control[1])")
print("vecS[1] = " + trD2[1] + ";")
print("if(control[2])")
print("vecS[2] = GSL_REAL(trH0);")
print("if(control[3])")
print("vecS[3] = GSL_REAL(trH0H0);")
print("if(control[4])")
print("vecS[4] = GSL_REAL(trH1);")
print("if(control[5])")
print("vecS[5] = GSL_REAL(trH1H1);")
print("}")

"""

#"""
#for simulation
#for D8+g6*D6+g4*D4+g2*D2

D8 = dirac_power(list1, 8)
D6 = dirac_power(list1, 6)
D4 = dirac_power(list1, 4)
D2 = dirac_power(list1, 2)

trD2 = masterchef(D2, init, init_mod, 0)
trD4 = masterchef(D4, init, init_mod, 0)
trD6 = masterchef(D6, init, init_mod, 0)
trD8 = masterchef(D8, init, init_mod, 1)


print("return G2*(" + trD2[0] + ") + G4*(" + trD4[0] + ") + G6*(" + trD6[0] + ") " + trD8[0] + ";")
print("}")
#"""

"""
#for D2
D2 = dirac_power(list1, 2)

trD2 = masterchef(D2, init, init_mod, 1)

print("return " + trD2[0] + ";")
print("}")
"""
