#!/usr/bin/python
# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: pyEPSim, a Python dynamical simulator for Electric Power Systems.
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 5.2
# DESCRIPTION: this program is a state estimator for an electric power system.
# The program is designed to dynamically simulate an Electric Power System provided the descriptive files of the system through integration of its Algebraic-Differential equations.
# BIBLIOGRAPHY:
# [1] 
# [2]
# [3]
# -------------------------------------------------

__author__ = "Alvaro Augusto Volpato and Luís Fernando Costa Alberto"
__license__ = "GPL"
__version__ = "5.2"
__maintainer__ = "Rob Knight"
__email__ = "alvaro.volpato@usp.br"
__status__ = "Development"

# -------------------------------------------------
# (1) IMPORTING LIBRARIES  {{{1
# -------------------------------------------------
# (1.1) Numpy libraries for standard operations  {{{2
import numpy as np
norm = np.linalg.norm
abs = np.absolute
cos = np.cos
sin = np.sin
array = np.array
real = np.real
imag = np.imag
conc = np.concatenate
transpose = np.transpose
inv = np.linalg.inv
sqrt = np.sqrt
pi = np.pi

# (1.2) Matplotlib {{{2
import matplotlib.pyplot as pyplt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
pyplt.rc('text', usetex=True)
pyplt.rc('text.latex', preamble=r'\usepackage{xfrac}')

# (1.2) Scipy for LU decomposition and ODE integration {{{2
import scipy
import scipy.linalg
import scipy.integrate
odeint = scipy.integrate.odeint
LU = scipy.linalg.lu

# (1.3) COPY library for temporary variable deep copy {{{2
from copy import deepcopy
copy = deepcopy

# (1.5) Time library: allows to compute execution time  {{{2
import time

# (1.6) Jacobian matrix functions {{{2
from libraries import matrixFunctions as mF

# (1.8) Classes library {{{2

import libraries.classes as cL		# classes.py contains the bus, branch, generator and case classes

import libraries.dynamicModels as dM	# dynamicModels.py contains the ODEs for loads and generators used

import libraries.arguments as pA		# arguments.py contains the available arguments passed to the program

import libraries.loadCase as lC		# loadCase.py contains the case-loading routine which reads the net file

# -------------------------------------------------

# -------------------------------------------------
# (3) DESCRIBRING SYSTEM: READING FILE {{{1
# -------------------------------------------------

args = pA.parseArguments()

netFile = args.net

case14 = lC.loadCase('case14.net')

nBus = case14.nBus

V, theta, res = pF.powerFlow(case14)

## -------------------------------------------------
## (11) STARTING NUMERICAL METHOD FOR DYNAMICAL SIMULATION {{{1
## -------------------------------------------------
#if dyn:
#	# Gathering simulation parameters from args
#	tFinal = args.time
#	pps = args.pointspersecond
#
#	print('\n' + '-'*50 + '\n DYNAMIC FAULT SIMULATION \n' + '-'*50)
#	print('\n --> Setting up dynamical simulation...')
#
#	# (11.1) Pre-fault state
#	YLoad = [(Pload[i] - 1j*Qload[i])/V[i]**2 for i in range(nBus)]	# YLoad is the equivalent conductance load matrix. Each load is modelled as a constant impedance
#	
#	# Permutating the Y and K matrixes. Explanation: the Y matrix built before did not consider which bars were generators and which are not. Nevertheless, to build the Y_RED matrix, the generator buses must be the first rows. This means that we must obtain a new Y matrix wherein the first n columns are the genertor buses.
#	# In order to do this, a permutation matrix PDyn is built. The building process works as follows:
#	PDyn = np.eye(nBus)
#		
#	genDataDyn = copy(genData)
#	busNameDyn = copy(busName)
#	for i in range(nGen):
#		for k in range(int(genDataDyn[i,0])):
#			if not mF.isGen(k,genDataDyn):
#				temp = copy(PDyn[int(genDataDyn[i,0]-1)])
#				PDyn[int(genDataDyn[i,0]-1)] = copy(PDyn[k])
#				PDyn[k] = temp
#				genDataDyn[i,0] = k+1
#				break
#
#	# Remember that we need to permute both columns and rows. To permute rows, it is needed to left-multiply Y by P (PY) and to columns swap, right-multiply (YP).
#	YDyn = PDyn@Y@PDyn
#	YLoadDyn = PDyn@YLoad
#	VDyn,thetaDyn = PDyn@V, PDyn@theta
#
#	Yred,C,D = mF.reduceGrid(YDyn,YLoadDyn,VDyn,genData)
#
#	
#	# Pre-fault mechanical power
#	pm = np.zeros(nGen)
#	for i in range(nGen):
#		pm[i] = (VDyn[i]**2)*real(Yred[i,i]) + sum( [ (C[i,j]*sin(thetaDyn[i] - thetaDyn[j]) + D[i,j]*cos(thetaDyn[i] - thetaDyn[j])) for j in range(nGen)])
#
#	fig = [ [] for i in range(nFault)]
#
#	for q in range(nFault):
#
#		lineAdmittance = Y[int(faultData[q,1]),int(faultData[q,2])] 
#
#		YFault = copy(Y)
#		YLoadFault = copy(YLoad)
#		VDynFault = copy(V)
#	
#		YFault[int(faultData[q,1]),int(faultData[q,2])] = 0
#		YFault[int(faultData[q,2]),int(faultData[q,1])] = 0
#
#
#		if int(faultData[q,3]) == 0:
#			#YLoadFault[int(faultData[q,1])] += lineAdmittance*1e3
#			YLoadFault[int(faultData[q,1])] = 0
#			VDynFault[int(faultData[q,1])] = 0
#			YLoadFault[int(faultData[q,2])] += lineAdmittance
#		elif int(faultData[q,3]) == 1:
#			YLoadFault[int(faultData[q,1])] += lineAdmittance
#			YLoadFault[int(faultData[q,2])] += lineAdmittance*1e3
#		else:
#			YLoadFault[int(faultData[q,1])] += lineAdmittance/faultData[q,3]
#			YLoadFault[int(faultData[q,2])] += lineAdmittance/(1-faultData[q,3])
#
#		YDynFault = PDyn@YFault@PDyn
#		YLoadDynFault = PDyn@YLoadFault
#		VDynFault,thetaDynFault = PDyn@VDynFault, PDyn@theta
#
#		YredFault,CFault,DFault = mF.reduceGrid(YDynFault,YLoadDynFault,VDynFault,genDataDyn)
#
#		print(' --> Starting numerical simulation of fault {0}...'.format(int(faultData[q,0] + 1)),end='')
#
#
#		tFault = np.linspace(0,faultData[q,4],int(round(pps*(faultData[q,4]))),endpoint=True)
#		y0 = np.zeros(2*nGen)
#
#		for i in range(nGen): y0[2*i] = thetaDyn[i]
#
#		solFault = odeint(mF.sm2,y0,tFault,args=(CFault,DFault,YredFault,VDyn,pm,genDataDyn))
#
#		# Calculating post-fault system
#
#		YPost = copy(YDyn)
#		YLoadPost = copy(YLoadDyn)
#
#		YPost[int(faultData[q,1]),int(faultData[q,2])] = 0
#		YPost[int(faultData[q,2]),int(faultData[q,1])] = 0
#
#		YDynPost = PDyn@YFault@PDyn
#		YLoadDynPost = PDyn@YLoadPost
#
#		YredPost,CPost,DPost = mF.reduceGrid(YDynPost,YLoadDynPost,VDyn,genDataDyn)
#
#		tPost = np.linspace(faultData[q,4],tFinal,int(round(pps*(tFinal - faultData[q,4]))),endpoint=True)
#		solPost = odeint(mF.sm2,solFault[-1,:],tPost,args=(CPost,DPost,YredPost,VDyn,pm,genDataDyn))
#
#		print('Done.')
#
#		# Plotting results
#
#		sol = conc((solFault,solPost),axis=0)
#		t = conc((tFault,tPost),axis=0)
#
#		fig[q] = pyplt.figure()
#		ax1 = fig[q].add_subplot(2,1,1)
#		ax2 = fig[q].add_subplot(2,1,2)
#
#		for i in range(nGen):
#			ax2.plot(t,sol[:, 2*i] - sol[:, 0])
#			ax1.plot(t,sol[:, 2*i+1]-sol[:, 1],label=busName[int(genData[i,0])-1])
#
#		ax1.axvline(x=faultData[q,4], color='k', lw=1, linestyle='--',label=r'$t_F$'.format(faultData[q,4]))
#		ax2.axvline(x=faultData[q,4], color='k', lw=1, linestyle='--',label=r'$t_F$'.format(faultData[q,4]))
#		ax1.legend(loc="lower right")
#
#		ax1.set_ylabel(r'Synchronized angular speed $\left(\sfrac{rad}{s}\right)$')
#		ax2.set_ylabel(r'Power angle $\left(rad\right)$')
#		ax1.grid(which='major',axis='both')
#		ax1.grid(which='major',axis='both')
#
#		ax1.grid(which='major',axis='both')
#		ax2.grid(which='major',axis='both')
#
#		fig[q].suptitle(r'Fault {0} results. Description: short-circuit at line connecting buses {1} and {2}, at position {3}, with open time {4} ms'.format(int(faultData[q,0])+1, int(faultData[q,1]+1), int(faultData[q,2]+1), faultData[q,3], 1000*faultData[q,4]))
#		ax2.set_xlabel('Time (s)')
#
#	print(' --> Plotting dynamical simulation results.')
#	pyplt.legend()
#
#	pyplt.show()
