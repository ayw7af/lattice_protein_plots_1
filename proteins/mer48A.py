#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Sequence of one 48-mer and inital configuration.
seq = "RQPDFYEQKDKTVERLRGGMGIATENYTNSASACILWHPFVDKLVSAL"


# Native list of the 48 mer.
# See reference:
# Shakhnovich, E., et. al. Science, 1996, 379, 96-98.
"""
(1, 4), (1, 40), (1, 42),
(2, 39),
(3, 6), (3, 8), 
(4, 9), (4, 15), 
(5, 16), (5, 18), (5, 20), (5, 40),
(6, 39), (6, 27),
(7, 18), (7, 28),
(9, 18),
(10, 15), (10, 17), 
(11, 14), 
(12, 17), (12, 33), 
(13, 16), (13, 44), (13, 48),
(14, 43), 
(15, 42), 
(16, 41), (16, 35),
(17, 34), 
(19, 28), (19, 30), (19, 34), 
(20, 27), (20, 35), (20, 37),  
(21, 24), (21, 26), (21, 30), 
(22, 31), (22, 47), (22, 35), 
(23, 36), (23, 46), 
(24, 37), 
(25, 38), 
(26, 29), 
(27, 38),
(31, 34), 
(32, 47), 
(33, 48), 
(35, 48), 
(36, 41), (36, 45), 
(37, 40), 
(41, 44), 
(45, 48)
"""

# Please note that the indexing of Python starts with 0.
native_list = [(0, 3), (0, 39), (0, 41), (1, 38), (2, 5), (2, 7), 
               (3, 8), (3, 14), (4, 39), (4, 15), (4, 19), (4, 17), 
               (5, 38), (5, 26), (6, 17), (6, 27), (8, 17), (9, 14), (9, 16),
               (10, 13), (11, 16), (11, 32), (12, 15), (12, 43), (12, 47), 
               (13, 42), (14, 41), (15, 40), (15, 34), (16, 33), 
               (18, 27), (18, 29), (18, 33), (19, 34), (19, 36), (19, 26), 
               (20, 23), (20, 29), (20, 25), (21, 30), (21, 46), (21, 34), 
               (22, 45), (22, 35), (23, 36), (24, 37), (25, 28), (26, 37), 
               (30, 33), (31, 46), (32, 47), (34, 47), (35, 44), (35, 40), 
               (36, 39), (40, 43), (44, 47)]