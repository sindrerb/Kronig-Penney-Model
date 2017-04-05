# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 15:01:28 2017

@author: sindrerb
"""
import numpy as np
import matplotlib.pyplot as plt

def FermiDirac(energy,fermilevel,temperature):
    return (np.exp((energy-fermilevel)/temperature)+1.0)**-1


basis  = np.loadtxt("./Kronig-Penney-Model/build/BASISFILE")
basisLength = len(basis[:,0])

""" SPLIT FILE INFORMATION """
states = np.loadtxt("./Kronig-Penney-Model/build/WAVEFILE")
kPoints = states[:,0:3]
energy = states[:,3]
waveFunction = states[:,4:]

""" CREATE BASIS MOMENTUM MATRIX """
G = np.zeros((basisLength,basisLength))
for i in range(0,basisLength):
    g = np.asarray(basis[i,:])
    G[i,i] = np.linalg.norm(g)

