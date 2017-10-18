import math
from numpy import matrix
from numpy import linalg
import numpy
import copy
import itertools




def generate_combinations(a, n):
    A = []
    for i in range(n):
        A.append(a)
    res1 = itertools.product(*A)
    res2 = [list(item) for item in res1]
    return res2


a = ["a", "b", "c", "d"]
print(len(generate_combinations(a, 4)))
