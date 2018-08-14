#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def abcd2pdb(abcd):
  ofilename = "latticechignolin.pdb"
  ofile = open(ofilename, "a+")
  x = []
  y = []
  z = []
  for i in range(len(abcd)+1):
    x.append(0.0)
    y.append(0.0)
    z.append(0.0)
  # print abcd
  for i in range(len(abcd)):
    if abcd[i]=="A":
      for j in range(i,len(abcd)):
        y[j+1] = y[j+1] - 1.4
    if abcd[i]=="B":
      for j in range(i,len(abcd)):
        x[j+1] = x[j+1] + 1.4
    if abcd[i]=="C":
      for j in range(i,len(abcd)):
        y[j+1] = y[j+1] + 1.4
    if abcd[i]=="D":
      for j in range(i,len(abcd)):
        x[j+1] = x[j+1] - 1.4
    if abcd[i]=="E":
      for j in range(i,len(abcd)):
        z[j+1] = z[j+1] + 1.4
    if abcd[i]=="F":
      for j in range(i,len(abcd)):
        z[j+1] = z[j+1] - 1.4
  #        P S V K M A
  ofile.write("ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f\n" % (x[0],y[0],z[0]))
  ofile.write("ATOM      2  CA  TYR     2    %8.3f%8.3f%8.3f\n" % (x[1],y[1],z[1]))
  ofile.write("ATOM      3  CA  ASP     3    %8.3f%8.3f%8.3f\n" % (x[2],y[2],z[2]))
  ofile.write("ATOM      4  CA  PRO     4    %8.3f%8.3f%8.3f\n" % (x[3],y[3],z[3]))
  ofile.write("ATOM      5  CA  GLU     5    %8.3f%8.3f%8.3f\n" % (x[4],y[4],z[4]))
  ofile.write("ATOM      6  CA  THR     6    %8.3f%8.3f%8.3f\n" % (x[5],y[5],z[5]))
  ofile.write("ATOM      7  CA  GLY     7    %8.3f%8.3f%8.3f\n" % (x[6],y[6],z[6]))
  ofile.write("ATOM      8  CA  THR     8    %8.3f%8.3f%8.3f\n" % (x[7],y[7],z[7]))
  ofile.write("ATOM      9  CA  TRP     9    %8.3f%8.3f%8.3f\n" % (x[8],y[8],z[8]))
  ofile.write("ATOM     10  CA  GLY    10    %8.3f%8.3f%8.3f\n" % (x[9],y[9],z[9]))
  ofile.write("END\n")
