#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import random

def changeit(conformation):
  #End moves
  
  if random.random()<0.5:
    r = random.randrange(6)
    if r==0:
      conformation = "A" + conformation[1:]
    if r==1:
      conformation = "B" + conformation[1:]
    if r==2:
      conformation = "C" + conformation[1:]
    if r==3:
      conformation = "D" + conformation[1:]
    if r==4:
      conformation = "E" + conformation[1:]
    if r==5:
      conformation = "F" + conformation[1:]
  if random.random()<0.5:
    r = random.randrange(6)
    if r==0:
      conformation = conformation[:-1] + "A"
    if r==1:
      conformation = conformation[:-1] + "B"
    if r==2:
      conformation = conformation[:-1] + "C"
    if r==3:
      conformation = conformation[:-1] + "D"
    if r==4:
      conformation = conformation[:-1] + "E"
    if r==5:
      conformation = conformation[:-1] + "F"
    #Corner moves
  for i in range(len(conformation)-1):
    part1 = conformation[i]
    part2 = conformation[i+1]
    if (part1=="A" or part1=="C") and (part2=="B" or part2=="D"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
    if (part1=="B" or part1=="D") and (part2=="A" or part2=="B"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
    if (part1=="B" or part1=="D") and (part2=="E" or part2=="F"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
    if (part1=="E" or part1=="F") and (part2=="B" or part2=="D"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
    if (part1=="A" or part1=="C") and (part2=="E" or part2=="F"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
    if (part1=="E" or part1=="F") and (part2=="A" or part2=="B"):
      if random.random()<0.5:
        conformation = conformation[:i]+part2+part1+conformation[i+2:]
  #Crankshaft moves
  for i in range(len(conformation)-2):
    part1 = conformation[i]
    part2 = conformation[i+1]
    part3 = conformation[i+2]
    if ((part1=="A" and part3=="C") or (part1=="C" and part3=="A")) and (part2=="B" or part2=="D"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
    if ((part1=="A" and part3=="C") or (part1=="C" and part3=="A")) and (part2=="E" or part2=="F"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
    if ((part1=="B" and part3=="D") or (part1=="D" and part3=="B")) and (part2=="A" or part2=="C"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
        conformation = conformation[:i]+part3+part2+part1+conformation[i+3:]
    if ((part1=="B" and part3=="D") or (part1=="D" and part3=="B")) and (part2=="E" or part2=="F"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
        conformation = conformation[:i]+part3+part2+part1+conformation[i+3:]
    if ((part1=="E" and part3=="F") or (part1=="F" and part3=="E")) and (part2=="A" or part2=="C"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"B"+part2+"D"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"D"+part2+"B"+conformation[i+3:]
        conformation = conformation[:i]+part3+part2+part1+conformation[i+3:]
    if ((part1=="E" and part3=="F") or (part1=="F" and part3=="E")) and (part2=="B" or part2=="D"):
      if random.random()<0.5:
        r = random.randrange(4)
        if r==0:
          conformation = conformation[:i]+"E"+part2+"F"+conformation[i+3:]
        if r==1:
          conformation = conformation[:i]+"F"+part2+"E"+conformation[i+3:]
        if r==2:
          conformation = conformation[:i]+"A"+part2+"C"+conformation[i+3:]
        if r==3:
          conformation = conformation[:i]+"C"+part2+"A"+conformation[i+3:]
        conformation = conformation[:i]+part3+part2+part1+conformation[i+3:]
  return conformation
  
  