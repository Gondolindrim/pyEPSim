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

def powerFlow(case,**kwargs):

	# absTol is the absolute tolerance used by the method
	if ('absTol' in kwargs): absTol = kwargs['absTol']
	else: absTol = 1e-6

	# deltaMax is the maximum dX at which the iterations are considered divergent. If dX > deltaMax the method aborts and returns success = False.
	if ('deltaMax' in kwargs): deltaMax = kwargs['deltaMax']
	else: deltaMax = 10

	# maxIter is the maximum number of iterations until the method is considered divergent.
	if ('maxIter' in kwargs): maxIter = kwargs['maxIter']
	else: maxIter = 10

	# Information print level.
	if ('verbose' in kwargs): verbose = kwargs['verbose']
	else: verbose = 0

	# Initial guess: flat start
	V = .9*np.ones(case.nBus)
	theta = np.zeros(case.nBus)

	# P and Q are the vectors of injected bus power
	P = np.diag(array([bus.pGen - bus.pLoad for bus in case.busData]))/case.Sb
	Q = np.diag(array([bus.qGen - bus.qLoad for bus in case.busData]))/case.Sb

	# isP, isQ and isV are the matrixes/array that flag P, Q and V measures
	isP = np.eye(case.nBus)
	isV = np.zeros(case.nBus)

	numP = int(isP.sum())
	numV = int(isV.sum())

	success = False

	itCount = 0

	# (5.6) Start time counte6
	tstart = time.time()

	# -------------------------------------------------
	# (10) STARTING NUMERICAL METHOD FOR POWER FLOW {{{1
	# -------------------------------------------------

	while(True):
	# (6) STARTING ITERATIONS

		# (6.1) Increasing iteration counter
		itCount += 1
		
		# (6.2) Printing iteration report
		if verbose > 0: print('\n ==> Iteration #{0:3.0f} '.format(itCount) + '-'*50)
		
		# (6.3) Calculating mF.Jacobian
		H = mF.Jac(V,theta,case.K,case.a,case.y,case.Y,case.bsh,isP,isV)
		H = np.delete(H,0,axis=0) # Removing slack bar angle derivatives for it is known
		Z = mF.Z(P,Q,isP,V,isV)
		h = mF.h(V,theta,case.K,case.a,case.y,case.Y,case.bsh,isP,isV)

		# (6.5) Calculating state update
		deltaSLC = np.delete(Z - h,0,axis=0)
		dX = inv(H) @ deltaSLC

		# (6.6) Updating V and theta
		theta[1:] += dX[0:case.nBus-1 ,0]
		V += dX[case.nBus-1:,0]
		
		# (6.8) Printing iteration results
		if verbose > 0: print(' --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}'.format(norm(dX),transpose(dX[case.nBus-1: ,0]),dX[0: case.nBus-1,0]))
		if verbose > 1: print('\n --> J = \n{0},\n\n --> r = Z - h = \n{2}*\n{1}'.format(H, (Z - h)/norm(Z-h), norm(Z-h)))

		# (6.9) Testing for iteration sucess or divergence
		if norm(dX) < absTol:
			success = True
			break	

		if norm(dX) > deltaMax: break
		if itCount > maxIter - 1: break

		# (6.10) Pausing for each iteration
		if verbose > 1: pause('\n --> Power flow method paused for next iteration. Press <ENTER> to continue.')


	# (6.11) Calculating elapsed time
	elapsed = time.time() - tstart

	return [V, theta, norm(dX), elapsed, itCount, success]
