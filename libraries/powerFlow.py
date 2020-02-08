# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: power flow method used by pyEPSim
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 16/08/2018
# VERSION: V1.0.0
# DESCRIPTION: this program is a power flow calculator for an electric power system. It receives a case object instance and keyword arguments which are explained below.
# -------------------------------------------------

import numpy as np
norm = np.linalg.norm
array = np.array
transpose = np.transpose
inv = np.linalg.inv

# Time library: allows to compute execution time 
import time

# Library functions
import libraries.classes as cL			# classes.py contains the bus, branch, generator and case classes
from libraries import matrixFunctions as mF	# Contains the matrix functions used to compute the Jacobian

# Tabulate for pretty tabular console output
from libraries.tabulate.tabulate import tabulate

from copy import deepcopy
copy = deepcopy

def powerFlow(case,**kwargs):

	# absTol is the absolute tolerance used by the method
	if ('absTol' in kwargs): absTol = kwargs['absTol']
	else: absTol = 1e-6

	# deltaMax is the maximum dX at which the iterations are considered divergent. If dX > deltaMax the method aborts and returns success = False.
	if ('deltaMax' in kwargs): deltaMax = kwargs['deltaMax']
	else: deltaMax = 100

	# maxIter is the maximum number of iterations until the method is considered divergent.
	if ('maxIter' in kwargs): maxIter = kwargs['maxIter']
	else: maxIter = 100

	# Information print level.
	# If verbose is 0, no output is given. This is the default level.
	# If 1, the results are output at the end of the method.
	# If 2, each iteration result is displayed and program is paused ultil a key is pressed.
	if ('verbose' in kwargs): verbose = kwargs['verbose']
	else: verbose = 0

	# Initial guess: flat start
	V = np.ones(case.nBus)	
	theta = np.zeros(case.nBus)

	# P and Q are the vectors of injected bus power
	P = np.diag(array([bus.pGen - bus.pLoad for bus in case.busData]))/case.Sb
	Q = np.diag(array([bus.qGen - bus.qLoad for bus in case.busData]))/case.Sb
	#print(P)
	# isP, isQ and isV are the matrixes/array that flag P, Q and V measures
	isP = np.eye(case.nBus)
	isQ = np.eye(case.nBus)
	isV = np.ones(case.nBus)
	isT = np.ones(case.nBus)
	
	for i in range(case.nBus):
		if case.busData[i].PVtype == 'VT':
			V[i] = case.busData[i].finalVoltage
			theta[i] = case.busData[i].finalAngle
			isV[i] = 0
			isT[i] = 0
			isP[i,i] = 0
			isQ[i,i] = 0
		if case.busData[i].PVtype == 'PV':
			V[i] = case.busData[i].finalVoltage
			isV[i] = 0
			isP[i,i] = 1
			isQ[i,i] = 0
		if case.busData[i].PVtype == 'PQ':
			isP[i,i] = 1
			isQ[i,i] = 1

	success = False

	itCount = 0
	verbose = 0
	# (5.6) Start time counte6
	tstart = time.time()
	# -------------------------------------------------
	# STARTING NUMERICAL METHOD FOR POWER FLOW
	# -------------------------------------------------
	if verbose > 0: print(' --> Beggining power flow method on case \'{0}\'...'.format(case.name))
	while(True):	# STARTING ITERATIONS

#		for i in range(case.nBus):
#			if case.busData[i].pLoad != 0 and case.busData[i].qLoad**2 != 0:
#				case.busData[i].gsh = -V[i]**2*case.busData[i].pLoad/(case.busData[i].pLoad**2 + case.busData[i].qLoad**2)
#				case.busData[i].bsh = -V[i]**2*case.busData[i].qLoad/(case.busData[i].pLoad**2 + case.busData[i].qLoad**2)
#		case.updateMatrixes()

		# Increasing iteration counter
		itCount += 1
		
		# Printing iteration report
		if verbose > 1: print('\n ==> Iteration #{0:3.0f} '.format(itCount) + '-'*50)
		
		# Calculating mF.Jacobian
		H = mF.Jac(V,theta,case.K,case.a,case.phi,case.y,case.Y,case.bsh,isP,isQ,isV,isT)	
		if verbose > 2: print(' >>>>> H = {}'.format(H))
		Z = mF.Z(P,Q,isP,isQ)
		h = mF.h(V,theta,case.K,case.a,case.phi,case.y,case.Y,case.bsh,isP,isQ,isV)

		# Calculating state update
		deltaSLC = Z - h

		dX = inv(H) @ deltaSLC
		
		# Updating V and theta
		i = 0
		for j in range(case.nBus):
			if isT[j]:
				theta[j] += dX[i]
				i+=1
		for j in range(case.nBus):
			if isV[j]:
				V[j] += dX[i]
				i+=1

		# Printing iteration results
		if verbose > 1: print(' --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}'.format(norm(dX),transpose(dX[case.nBus-1: ,0]),dX[0: case.nBus-1,0]))

		# Testing for iteration sucess or divergence
		if norm(dX) < absTol:
			success = True
			break	

		if norm(dX) > deltaMax: break
		if itCount > maxIter - 1: break

		# Pausing for each iteration
		if verbose > 1: input('\n --> Power flow method paused for next iteration. Press <ENTER> to continue.')

	# Calculating elapsed time
	elapsed = time.time() - tstart

	# Printing results
	if verbose > 0:

		tabHeaders = ['Bus name','Number','Voltage (p.u.)','Angle (rad)']

		tabRows = []
		for bus in case.busData:
			tabRows.append([bus.name, bus.number, V[bus.number], theta[bus.number]])
			
		resultsTable = tabulate(tabRows,headers=tabHeaders,tablefmt='psql')

		# Printing convergence results
		if success:
			print('\n >> Power flow method successful with {0} total iterations and step increment norm |dX| = {1}. Results:'.format(itCount,norm(dX)))
			print(resultsTable)
		#	if verbose > 1: print('\n --> Residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,r))

		else: print(' --> Power flow method not successful at iteration {0} with step increment |dX| = {1}.'.format(itCount,norm(dX)))

	return [V, theta, norm(dX), elapsed, itCount, success]
