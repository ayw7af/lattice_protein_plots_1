#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Begins changeit.py.
#
# Procedure to evolve conformations of a lattice protein.
# Please note: the conformations are intependent to sequence.
#
# In the future OOP, this will become a member function of the ``conf'' class.
# =============================================================================

import random
from copy import deepcopy
import numpy as np

def isselfavoiding(x, y, z, j, minus, plus):
    """
    Function that performs self-avoiding checks of residue index j.

    Returns True if the input conformation at index j satisfies self-avoiding,\
        returns False otherwise.
    """
    N = len(x)

    orientations = {"A": np.array([0, -1, 0]), "B": np.array([1, 0, 0]), \
                    "C": np.array([0, 1, 0]), "D": np.array([-1, 0, 0]), \
                    "E": np.array([0, 0, 1]), "F": np.array([0, 0, -1])}

    temp_x = deepcopy(x)
    temp_y = deepcopy(y)
    temp_z = deepcopy(z)

    change_array = orientations[plus] - orientations[minus]

    temp_x[j] += change_array[0]
    temp_y[j] += change_array[1]
    temp_z[j] += change_array[2]

    # Self-avoiding checks.
    for i in range(0, j - 1):
#        print(i)
        if (abs(temp_x[i]-temp_x[j])**2 + abs(temp_y[i]-temp_y[j])**2 + abs(temp_z[i]-temp_z[j])**2)==0:
#            print(False)
            return False
    for i in range(j + 2, N):
#        print(i)
        if (abs(temp_x[i]-temp_x[j])**2 + abs(temp_y[i]-temp_y[j])**2 + abs(temp_z[i]-temp_z[j])**2)==0:
#            print(False)
            return False

    return True

def isselfavoiding0(conformation):
    """
    ***DEPRECIATED METHOD***

    Function that performs self-avoiding checks of residue index k.

    Returns True if the input conformation at index k satisfies self-avoiding,\
        returns False otherwise.
    """
    x = []
    y = []
    z = []
    for i in range(len(conformation)+1):
        x.append(0)
        y.append(0)
        z.append(0)
    for i in range(len(conformation)):
        if conformation[i]=="A":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] - 1
        if conformation[i]=="B":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] + 1
        if conformation[i]=="C":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] + 1
        if conformation[i]=="D":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] - 1
        if conformation[i]=="E":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] + 1
        if conformation[i]=="F":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] - 1

    # Self-avoiding checks.
    for i in range(len(conformation)-1):
        for j in range(i+2, len(conformation)+1):
            if (abs(x[i]-x[j])**2 + abs(y[i]-y[j])**2 + abs(z[i]-z[j])**2)==0:
                return False

    return True

def changeit(conformation):
    """
    Procedure that performs conformational change.
    Only *** one *** conformational change will be produced in each \
        Monte Carlo step.

    Move set: Ends move, corners flip, and crankshaft change.
    """
    direction = ("A", "B", "C", "D", "E", "F")

    x = []
    y = []
    z = []
    for i in range(len(conformation)+1):
        x.append(0)
        y.append(0)
        z.append(0)
    for i in range(len(conformation)):
        if conformation[i]=="A":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] - 1
        if conformation[i]=="B":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] + 1
        if conformation[i]=="C":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] + 1
        if conformation[i]=="D":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] - 1
        if conformation[i]=="E":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] + 1
        if conformation[i]=="F":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] - 1

    # Count the number of allowed moves.
    # N_max(1) = N + 6.
    # N_max(2) = 3 * int((N - 2)/2).
    # N_max = N_max(1) + N_max(2) = 3 * int(N/2) + N + 3.
    N = len(conformation) + 1
    num_max = (N//2) * 3 + N + 3

    # Count the number of ends move.
    head_list = []
    for i in direction:
        if i == conformation[0]:
            continue
#        temp_conf = i + conformation[1:]
        if isselfavoiding(x, y, z, 0, i, conformation[0]) == True:
            head_list.append(i)
    num_head = len(head_list)

    tail_list = []
    for i in direction:
        if i == conformation[-1]:
            continue
#        temp_conf = conformation[:-1] + i
        if isselfavoiding(x, y, z, N - 1, conformation[-1], i) == True:
            tail_list.append(i)
    num_tail = len(tail_list)

    # Count the number of corners flip.
    # Corners_list stores the first index where it could perform a corner flip.
    corners_list = []
    for i in range(N - 2):
        part1 = conformation[i]
        part2 = conformation[i+1]
        if ((part1=="A" or part1=="C") and (part2=="B" or part2=="D")) or \
           ((part1=="B" or part1=="D") and (part2=="A" or part2=="C")) or \
           ((part1=="B" or part1=="D") and (part2=="E" or part2=="F")) or \
           ((part1=="E" or part1=="F") and (part2=="B" or part2=="D")) or \
           ((part1=="A" or part1=="C") and (part2=="E" or part2=="F")) or \
           ((part1=="E" or part1=="F") and (part2=="A" or part2=="C")):
#               temp_conf = conformation[:i] + part2 + part1 + conformation[i+2:]
               if isselfavoiding(x, y, z, i+1, part1, part2) == True:
                   corners_list.append(i)
        num_corners = len(corners_list)

    # Count the number of crankshaft change.
    crankshaft_list = []
    for i in range(N - 3):
        part1 = conformation[i]
        part2 = conformation[i+1]
        part3 = conformation[i+2]

        if ((part1=="A" and part3=="C") or (part1=="C" and part3=="A")) \
            and (part2=="B" or part2=="D"):
                if part1 =="C":
#                    temp_conf = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "C", "A") == True and \
                       isselfavoiding(x, y, z, i+2, "C", "A") == True:
                           crankshaft_list.append((i, "A", "C"))
                else:
#                    temp_conf = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "A", "C") == True and \
                       isselfavoiding(x, y, z, i+2, "A", "C") == True:
                           crankshaft_list.append((i, "C", "A"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "E") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "E") == True:
                    crankshaft_list.append((i, "E", "F"))
#                temp_conf = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1 ,part1, "F") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "F") == True:
                    crankshaft_list.append((i, "F", "E"))

        elif ((part1=="A" and part3=="C") or (part1=="C" and part3=="A")) \
            and (part2=="E" or part2=="F"):
                if part1 =="C":
#                    temp_conf = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "C", "A") == True and \
                       isselfavoiding(x, y, z, i+2, "C", "A") == True:
                           crankshaft_list.append((i, "A", "C"))
                else:
#                    temp_conf = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "A", "C") == True and \
                       isselfavoiding(x, y, z, i+2, "A", "C") == True:
                           crankshaft_list.append((i, "C", "A"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "B") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "B") == True:
                       crankshaft_list.append((i, "B", "D"))
#                temp_conf = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "D") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "D") == True:
                       crankshaft_list.append((i, "D", "B"))

        elif ((part1=="B" and part3=="D") or (part1=="D" and part3=="B")) \
            and (part2=="A" or part2=="C"):
                if part1 =="D":
#                    temp_conf = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "D", "B") == True and \
                       isselfavoiding(x, y, z, i+2, "D", "B") == True:
                           crankshaft_list.append((i, "B", "D"))
                else:
#                    temp_conf = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "B", "D") == True and \
                       isselfavoiding(x, y, z, i+2, "B", "D") == True:
                           crankshaft_list.append((i, "D", "B"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "E") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "E") == True:
                    crankshaft_list.append((i, "E", "F"))
#                temp_conf = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "F") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "F") == True:
                    crankshaft_list.append((i, "F", "E"))

        elif ((part1=="B" and part3=="D") or (part1=="D" and part3=="B")) \
            and (part2=="E" or part2=="F"):
                if part1 =="D":
#                    temp_conf = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "D", "B") == True and \
                       isselfavoiding(x, y, z, i+2, "D", "B") == True:
                           crankshaft_list.append((i, "B", "D"))
                else:
#                    temp_conf = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "B", "D") == True and \
                       isselfavoiding(x, y, z, i+2, "B", "D") == True:
                           crankshaft_list.append((i, "D", "B"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "A") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "A") == True:
                    crankshaft_list.append((i, "A", "C"))
#                temp_conf = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "C") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "C") == True:
                    crankshaft_list.append((i, "C", "A"))

        elif ((part1=="E" and part3=="F") or (part1=="F" and part3=="E")) \
            and (part2=="A" or part2=="C"):
                if part1 =="F":
#                    temp_conf = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "F", "E") == True and \
                       isselfavoiding(x, y, z, i+2, "F", "E") == True:
                           crankshaft_list.append((i, "E", "F"))
                else:
#                    temp_conf = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "E", "F") == True and \
                       isselfavoiding(x, y, z, i+2, "E", "F") == True:
                           crankshaft_list.append((i, "F", "E"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "B") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "B") == True:
                       crankshaft_list.append((i, "B", "D"))
#                temp_conf = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "D") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "D") == True:
                       crankshaft_list.append((i, "D", "B"))

        if ((part1=="E" and part3=="F") or (part1=="F" and part3=="E")) \
            and (part2=="B" or part2=="D"):
                if part1 =="F":
#                    temp_conf = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "F", "E") == True and \
                       isselfavoiding(x, y, z, i+2, "F", "E") == True:
                           crankshaft_list.append((i, "E", "F"))
                else:
#                    temp_conf = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
                    if isselfavoiding(x, y, z, i+1, "E", "F") == True and \
                       isselfavoiding(x, y, z, i+2, "E", "F") == True:
                           crankshaft_list.append((i, "F", "E"))
#                if isselfavoiding(x, y, z, temp_conf) == True:
#                    crankshaft_list.append((i, temp_conf[i], temp_conf[i+2]))
#                temp_conf = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "A") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "A") == True:
                    crankshaft_list.append((i, "A", "C"))
#                temp_conf = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
                if isselfavoiding(x, y, z, i+1, part1, "C") == True and \
                   isselfavoiding(x, y, z, i+2, part1, "C") == True:
                    crankshaft_list.append((i, "C", "A"))

    num_crankshaft = len(crankshaft_list)

    # Perform move.
    """
    Creates a random number that within the range of number of all moves.
    if r < num_ends, perform corresponding ends move.
    if num_ends < r < num_corners, perform corner flip.
    otherwise, perform crankshaft change.
    """
    r = random.randrange(num_max)

    # Ends move.
    if r < num_head:
        conformation = head_list[r] + conformation[1:]
        return conformation
    elif r < num_head + num_tail:
        conformation = conformation[:-1] + tail_list[r - num_head]
        return conformation
    elif r < num_head + num_tail + num_corners:
        # Locate the ID of atom that performs corner flip.
        i = corners_list[r - num_head - num_tail]
        part1 = conformation[i]
        part2 = conformation[i+1]
        conformation = conformation[:i] + part2 + part1 + conformation[i+2:]
        return conformation
    elif r < num_head + num_tail + num_corners + num_crankshaft:
       # Locate the ID of atom that performs crankshaft change.
        (i, part1, part3) = crankshaft_list[r - num_head - num_tail - num_corners]
        conformation = conformation[:i] + part1 + conformation[i+1] + \
                       part3 + conformation[i+3:]
        return conformation

    return conformation