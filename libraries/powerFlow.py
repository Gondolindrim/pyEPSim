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
	else: deltaMax = 10

	# maxIter is the maximum number of iterations until the method is considered divergent.
	if ('maxIter' in kwargs): maxIter = kwargs['maxIter']
	else: maxIter = 10

	# Information print level.
	# If verbose is 0, no output is given. This is the default level.
	# If 1, the results are output at the end of the method.
	# If 2, each iteration result is displayed and program is paused ultil a key is pressed.
	if ('verbose' in kwargs): verbose = kwargs['verbose']
	else: verbose = 0

	# fixedAngleList is the list of buses that have a fixed angle and must not be changed during the powerFlow method.
	for bus in case.busData: 
		if bus.PVtype == 'VT':
			 slackBusNumber = bus.number

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
	# STARTING NUMERICAL METHOD FOR POWER FLOW
	# -------------------------------------------------
	if verbose > 0: print(' --> Beggining power flow method on case \'{0}\'...'.format(case.name))
	while(True):	# STARTING ITERATIONS

		# Increasing iteration counter
		itCount += 1
		
		# Printing iteration report
		if verbose > 1: print('\n ==> Iteration #{0:3.0f} '.format(itCount) + '-'*50)
		
		# Calculating mF.Jacobian
		H = mF.Jac(V,theta,case.K,case.a,case.y,case.Y,case.bsh,isP,isV)
		H = np.delete(H,0,axis=0) # Removing slack bar angle derivatives for it is known
		Z = mF.Z(P,Q,isP,V,isV)
		h = mF.h(V,theta,case.K,case.a,case.y,case.Y,case.bsh,isP,isV)

		# Calculating state update
		deltaSLC = np.delete(Z - h,0,axis=0)
		dX = inv(H) @ deltaSLC
	
		# Updating V and theta
		theta[1:] += dX[0:case.nBus-1 ,0]
		V += dX[case.nBus-1:,0]
		
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
