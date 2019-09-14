########## maxcut.py ##########

import localsolver
import sys

# if len(sys.argv) < 2:
#     print ("Usage: python maxcut.py inputFile [outputFile] [timeLimit]")
#     sys.exit(1)


def read_integers(filename):
    with open(filename) as f:
        return [int(elem) for elem in f.read().split()]


with localsolver.LocalSolver() as ls:

    #
    # Reads instance data
    #

    file_it = iter(read_integers(sys.argv[1]))
    # Number of vertices
    n = next(file_it)
    # Number of edges
    m = next(file_it)

    # Origin of each edge
    origin = [None]*m
    # Destination of each edge
    dest = [None]*m
    # Weight of each edge
    w = [None]*m

    for e in range(m):
        origin[e] = next(file_it)
        dest[e] = next(file_it)
        w[e] = next(file_it)

    #
    # Declares the optimization model
    #
    model = ls.model

    # Decision variables x[i]
    # Is true if vertex x[i] is on the right side of the cut and false if it is on the left side of the cut
    x = [model.bool() for i in range(n)]

    # incut[e] is true if its endpoints are in different class of the partition
    incut = [None]*m
    for e in range(m):
        incut[e] = model.neq(x[origin[e] - 1], x[dest[e] - 1])

    # Size of the cut
    cut_weight = model.sum(w[e]*incut[e] for e in range(m))
    model.maximize(cut_weight)

    model.close()

    #
    # Param
    #
    # if len(sys.argv) >= 4: ls.param.time_limit = int(x3)
    # else: ls.param.time_limit = 10
    ls.param.time_limit = 5
    ls.solve()

    #
    # Writes the solution in a file following the following format:
    #  - objective value
    #  - each line contains a vertex number and its subset (1 for S, 0 for V-S)
    #
    # if len(sys.argv) >= 3:
    with open(sys.argv[2], 'w') as f:
        f.write("%d %d\n" % (cut_weight.value,1))
        # Note: in the instances the indices start at 1
        for i in range(n):
            f.write("%d %d\n" % (i+1, x[i].value))
