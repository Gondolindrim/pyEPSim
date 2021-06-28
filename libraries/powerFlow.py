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

def power_flow(case,**kwargs): #{{{1
	print(case.isP)
	# absTol is the absolute tolerance used by the method
	if ('absTol' in kwargs): absTol = kwargs['absTol']
	else: absTol = 1e-6

	# deltaMax is the maximum dX at which the iterations are considered divergent. If dX > deltaMax the method aborts and returns success = False.
	if ('deltaMax' in kwargs): deltaMax = kwargs['deltaMax']
	else: deltaMax = 1e20

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
	nbus = case.get_nbus()
	V = np.ones(nbus)	
	theta = np.zeros(nbus)
	for n in range(nbus): 
		if case.bus_data[n].bus_type in ['VT', 'PV']: V[n] = case.bus_data[n].final_voltage
		if case.bus_data[n].bus_type == 'VT': theta[n] = case.bus_data[n].final_angle
	
	# P and Q are the vectors of injected bus power. They are calculated as the difference between the injected power on the bus and the load. This difference does not account for 
	#	shunt impedance loads, only for power (MW and MVAr) loads. The impendance loads are dealth with through the impedance matrices.
	P = np.zeros((nbus, nbus))
	Q = np.zeros((nbus, nbus))
	for i in range(nbus):
		bus = case.bus_data[i]
		gen = cL.is_gen(bus.name, case.gen_data)
		if gen:
			P[i,i] = gen.p_gen
			Q[i,i] = gen.q_gen

		P[i,i] = P[i,i] - bus.p_load
		Q[i,i] = Q[i,i] - bus.q_load

	P, Q = P/case.Sb, Q/case.Sb

	success = False
	iteration_count = 0
	verbose = 0

	# -------------------------------------------------
	# STARTING NUMERICAL METHOD FOR POWER FLOW
	# -------------------------------------------------
	if verbose > 0: print(' --> Beggining power flow method on case \'{0}\'...'.format(case.name))
	while(True):	# STARTING ITERATIONS
		# Increasing iteration counter
		iteration_count += 1
		
		# Printing iteration report
		if verbose > 1: print('\n ==> Iteration #{0:3.0f} '.format(iteration_count) + '-'*50)
		
		# Calculating mF.Jacobian
		H = mF.Jac(V,theta,case)	
		if verbose > 2: print(' >>>>> H = {}'.format(H))
		Z = mF.Z(P,Q,case)
		h = mF.h(V,theta,case)
		# Calculating state update
		deltaSLC = Z - h

		dX = inv(H) @ deltaSLC
		
		# Updating V and theta
		i = 0
		for j in range(nbus):
			if case.isT[j]:
				theta[j] += dX[i]
				i+=1
		for j in range(nbus):
			if case.isV[j]:
				V[j] += dX[i]
				i+=1
			
		# Printing iteration results
		if verbose > 1: print(' --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}'.format(norm(dX),transpose(dX[case.nBus-1: ,0]),dX[0: case.nBus-1,0]))

		# Testing for iteration sucess or divergence
		if norm(dX) < absTol:
			success = True
			break	

		if norm(dX) > deltaMax: break
		if iteration_count > maxIter - 1: break

		# Pausing for each iteration
		if verbose > 1: input('\n --> Power flow method paused for next iteration. Press <ENTER> to continue.')

	# Printing results
	if verbose > 0:

		tab_headers = ['Bus name','Number','Voltage (p.u.)','Angle (rad)']

		tab_rows = []
		for bus in case.bus_data:
			tab_rows.append([bus.name, bus.number, V[bus.number], theta[bus.number]])
			
		results_table = tabulate(tab_rows, headers = tab_headers, tablefmt='psql')

		# Printing convergence results
		if success:
			print('\n >> Power flow method successful with {0} total iterations and step increment norm |dX| = {1}. Results:'.format(iteration_count,norm(dX)))
			print(results_table)
		#	if verbose > 1: print('\n --> Residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,r))

		else: print(' --> Power flow method not successful at iteration {0} with step increment |dX| = {1}.'.format(iteration_count,norm(dX)))

	return [V, theta, norm(dX), iteration_count, success]
