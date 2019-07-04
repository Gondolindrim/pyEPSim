#!/usr/bin/python
# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# AUTHOR: Álvaro Augusto Volpato
# DATE: 04/07/2018
# VERSION: 5.0
# DESCRIPTION: this program is a state estimator for an electric power system.
# The program first tests for observability, restores it from available pseudo-measures if needed, and pro-
# ceeds to find critical measurements. Next the program estimates the network's states through a Newton-
# -Raphson method proposed in MONTICELLI(1995), and tests for normalized residue.
# -------------------------------------------------

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

import matplotlib.pyplot as pyplt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
pyplt.rc('text', usetex=True)
pyplt.rc('text.latex', preamble=r'\usepackage{xfrac}')


# (1.2) SCIPY FOR LU DECOMPOSITION 
import scipy
import scipy.linalg
import scipy.integrate
odeint = scipy.integrate.odeint
LU = scipy.linalg.lu

# (1.3) COPY library for temporary variable deep copy 
from copy import deepcopy
copy = deepcopy

# (1.4) OS and SYS libraries for defining clear terminal command 
import os
import sys

# (1.5) Time library: allows to count execution time 
import time

# (1.6) Argument Parser to pass arguments to script  2
import argparse
parser = argparse.ArgumentParser(
	prog = 'Power flow calculator',
	description = 'Power flow calculator that takes a system data in the IEEE CDF and a measurements file in the LACO CDF.',
	)

# (1.7) mF.Jacobian matrix functions
import matrixFunctions as mF

#import linecache as lc

# -------------------------------------------------
# (2) CUSTOM FUNCTIONS AND ARGUMENTS{{{1
# -------------------------------------------------
# (2.1) Set printout options: this allows the program printout to span the whole terminal 
np.set_printoptions(suppress=False,linewidth=999,threshold=999)

# (2.2) Clear  function: clears the terminal. Similar to MATLAB's cls and bash's cls
def clear(): os.system('cls' if os.name == 'nt' else 'clear')

# (2.3) "pause" function: halts program execution and asks for enter 
def pause(string):
    programPause = input(string)

# (2.4) readtabline reads a line in a given file f then arranges it in an array which vectors are the strings
#	in the line that are separated by the tab character '\t'. It also removes the new line '\n' character
#	from the last slot in the line.
def readtabline(f): return f.readline().strip().split('\t')

# (2.5) Arguments and options 
#	Creates arguments:
#	--net [NETWORK FILE] is the target system network file. This string is stored in a 'networkFile' variable
parser.add_argument(	"net",
			type = str,
			metavar = 'NETFILE',
			help='Target network file')

# VERBOSE option
parser.add_argument(	'--verbose', '-v',
			type = int,
			choices = [0,1,2],
			default = [1],
			nargs = 1,
			help = 'Verbose level. 0 for no output, 1 for critical output, 2 for detailed output. Default is 1.')

# OBS option
parser.add_argument(	'--obs',
			default = False,
			action = 'store_true',
			help = 'Performs observability check. If the system is deemed not observable, the program tries to add the supplied pseudo-measures from the measurements file. Default value is false.')

# CRIT option
parser.add_argument(	'--crit',
			default = False,
			action = 'store_true',
			help = 'Performs critical measurements check. Default is false.')

# CLEAR option
parser.add_argument(	'--cls',
			default = False,
			action = 'store_true',
			help = 'Clears the terminal before executing the program.')

# TOL option
parser.add_argument(	'--tol',
			type=float,
			default = 1e-6,
			help = 'Newton-Raphson numerical method tolerance. Default is 1e-6.')

# MAXITER option
parser.add_argument(	'--maxiter',
			type = int,
			default = 10,
			help = 'Newton-Raphson numerical method maximum iterations. Default is 10.')

# DELTAMAX OPTION
parser.add_argument(	'--deltamax',
			type = int,
			default = 10,
			help = 'Newton-Raphson numerical method maximum  step magnitude to consider as divergence. Default is 10.')

# POWER FLOW OPTION
parser.add_argument(	'--powerflow', '-pf',
			default = False,
			action = 'store_true',
			help = 'Perform Power Flow calculation of the target system.')

# MEASUREMENTS FILE
parser.add_argument(	'--measfile', '-meas',
			type = str,
			metavar = 'MEASFILE',
			help = 'Measurements file for state estimation.')

# STATE ESTIMATION
parser.add_argument(	'--estimatestate', '-se',
			default = False,
			action = 'store_true',
			help = 'Perform State Estimation of the Target Sysem. A measurement file needs to be provided.')

# DYNAMICAL FAULT SIMULATION
parser.add_argument(	'--dynamicalfault', '-dyn',
			action = 'store_true',
			default = False,
			help = 'Perform dynamical simulation of the faults in the system. These faults should be listed in the MEAS file.')

# FINAL SIMULATION TIME
parser.add_argument(	'--time', '-t',
			type = float,
			default = 1,
			help = 'Dynamical simulation time length in seconds. Default is 1.')

# POINTS PER SECOND
parser.add_argument(	'--pointspersecond', '-pps',
			type = int,
			default = 1000,
			help = 'Number of points per second of the simulation plots. Default is 1000.')

# Gathering arguments and options
args = parser.parse_args()
cls = args.cls
if cls:
	clear()
	clear()

verbose = args.verbose[0]
measurementsFile = args.measfile
networkFile = args.net
obs = args.obs
crit = args.crit
powerflow = args.powerflow
estimatestate = args.estimatestate
dyn = args.dynamicalfault

if obs and not measurementsFile:
	sys.exit('--> Fatal error: observability analysis was requested but no measurements file was supplied.')
if crit and not measurementsFile:
	sys.exit('--> Fatal error: observability analysis was requested but no measurements file was supplied.')

# -------------------------------------------------
# (3) DESCRIBRING SYSTEM: READING FILE {{{1
# -------------------------------------------------

# (3.1) lineNumber returns the line number of the current pointer in a given file
def lineNumber(f):
	i = 0
	pointer = f.tell()	# Captures the pointer of the present line
	f.seek(0)		# Goes to the beggining of the file
	while f.tell() != pointer:
		# This loop sweeps the lines in the file and compares it to the stored pointer If there is a match, break and return the counter
	#	f.readline
	#	pointerTemp = f.tell()
	#	if pointerTemp == pointer: break
		i += 1
		f.readline()

	return i

print('\n' + '-'*50 + '\n READING FILES \n' + '-'*50)

# (3.1) READING NETWORK FILE AND BUILDING TARGET SYSTEM DATA {{{2
with open(networkFile,'r') as fileData:

	breakFlag = True

	# (3.1.2) Acquiring bus data parameters
	print('\n --> Beginning bus data reading... ',end='' )

	# Searching for bus data start card
	while breakFlag==True:
		line = readtabline(fileData)
		if line[0] == 'MVA base:':
			Sb = float(line[1])	# Extracting power base value
			breakFlag = False

	# Searching for bus data start card
	breakFlag = True
	while breakFlag==True:
		line = readtabline(fileData)
		if line[0] == 'BUS DATA FOLLOWS':
			breakFlag = False

	fileData.readline()	# Skip the '---' line		

	# Getting number of buses nBus
	startDataCardLineNumber = lineNumber(fileData)	# Store the line number of the start data card
	startDataCardLinePointer = fileData.tell()	# Store the pointer to this line as it will be needed later
	breakFlag = True
	while breakFlag == True:
		line = readtabline(fileData)
		if line[0] == '-999':			# Read lines until the -999 string is detected
			endDataCardLineNumber = lineNumber(fileData)	# Store the line number of the end data card
			nBus = endDataCardLineNumber - startDataCardLineNumber - 1	# Compute number of buses
			breakFlag  = False

	# Pre-allocating the matrixes
	bsh = np.zeros((nBus,nBus))
	gsh = np.zeros((nBus,nBus))
	Pload = np.zeros(nBus)
	Qload = np.zeros(nBus)
	Pger = np.zeros(nBus)
	Qger = np.zeros(nBus)
	busName = [[] for i in range(nBus)]
			
	fileData.seek(startDataCardLinePointer)	# Go back to the start data card line

	for i in range(nBus):
		line = readtabline(fileData)
		busName[i] = str(line[1])
		Pload[i] = float(line[7])/Sb
		Qload[i] = float(line[8])/Sb
		Pger[i] = float(line[9])/Sb
		Qger[i] = float(line[10])/Sb
		gsh[i,i] = float(line[15])
		bsh[i,i] = float(line[16])


	# (3.1.3) Searching for branch data start card
	print(' Done.\n --> Beginning branch data reading... ',end='' )

	breakFlag = True
	while breakFlag==True:
		line = readtabline(fileData)
		if line[0] == 'BRANCH DATA FOLLOWS':
			breakFlag = False

	fileData.readline()	# Skip the '---' line	

	# Getting number of branches nBranch
	startDataCardLineNumber = lineNumber(fileData)	# Store the line number of the start data card
	startDataCardLinePointer = fileData.tell()	# Store the pointer to this line as it will be needed later
	breakFlag = True
	while breakFlag == True:
		line = readtabline(fileData)
		if line[0] == '-999':			# Read lines until the -999 string is detected
			endDataCardLineNumber = lineNumber(fileData)	# Store the line number of the end data card
			nBranch = endDataCardLineNumber - startDataCardLineNumber - 1	# Compute number of buses
			breakFlag  = False

	# (3.1.4) Acquiring branch data parameters
	
	# Initializing resistance and reactance matrixes
	r = np.zeros((nBus,nBus))
	x = np.zeros((nBus,nBus))

	# Initializing tap levels matrix
	a = np.ones((nBus,nBus))
	
	# Initializing connection matrix
	K = [ [] for i in range(nBus)]

	fileData.seek(startDataCardLinePointer)	# Go back to the start data card line

	for i in range(nBranch):
		line = readtabline(fileData)
		fromBus = int(line[0])
		toBus = int(line[1])

		K[fromBus-1] += [toBus - 1]
		K[toBus-1] += [fromBus - 1]

		r[fromBus-1,toBus-1] = float(line[6])
		r[toBus-1,fromBus-1] = float(line[6])

		x[fromBus-1,toBus-1] = float(line[7])
		x[toBus-1,fromBus-1] = float(line[7])

		bsh[fromBus-1,toBus-1] = float(line[8])	
		bsh[toBus-1,fromBus-1] = float(line[8])

		if float(line[14]) != 0:
			a[fromBus-1,toBus-1] = 1/float(line[14])

	# (3.1.5) Building Y matrix
	print(' Done.\n --> Building Y matrix... ',end='')
	
	z = array([[r[m,n] + 1j*x[m,n] for m in range(nBus)] for n in range(nBus)])
	y = array([[1/z[m,n] if z[m,n] != 0 else 0 for m in range(nBus)] for n in range(nBus)])

	g = real(copy(y))
	b = imag(copy(y))

	bsh = array([[ bsh[m,n]/2 if m!=n else bsh[m,n] for n in range(nBus)] for m in range(nBus)])

	Y  = -a*transpose(a)*y

	for k in range(nBus): Y[k,k] = 1j*bsh[k,k] + sum([a[k,m]**2*y[k,m] + 1j*bsh[k,m] for m in K[k]])

	# (3.1.3) Searching for generator data start card
	print(' Done.\n --> Beginning generator data reading... ',end='' )

	breakFlag = True
	while breakFlag==True:
		line = readtabline(fileData)
		if line[0] == 'GENERATOR DATA FOLLOWS':
			breakFlag = False

	fileData.readline()	# Skip the '---' line	

	# Getting number of branches nBranch
	startDataCardLineNumber = lineNumber(fileData)	# Store the line number of the start data card
	startDataCardLinePointer = fileData.tell()	# Store the pointer to this line as it will be needed later
	breakFlag = True
	while breakFlag == True:
		line = readtabline(fileData)
		if line[0] == '-999':			# Read lines until the -999 string is detected
			endDataCardLineNumber = lineNumber(fileData)	# Store the line number of the end data card
			nGen = endDataCardLineNumber - startDataCardLineNumber - 1	# Compute number of buses
			breakFlag  = False

	genData = np.zeros((nGen,16))

	fileData.seek(startDataCardLinePointer)	# Go back to the start data card line

	for i in range(nGen):
		line = readtabline(fileData)
		for k in range(16):
			genData[i,k] = line[k]

	genData[:,2] /= Sb	
	# (3.1.4) Searching for fault data start card
	print(' Done.\n --> Beginning fault data reading... ',end='' )

	breakFlag = True
	while breakFlag==True:
		line = readtabline(fileData)
		if line[0] == 'FAULT DATA FOLLOWS':
			breakFlag = False

	fileData.readline()	# Skip the '---' line	

	# Getting number of branches nBranch
	startDataCardLineNumber = lineNumber(fileData)	# Store the line number of the start data card
	startDataCardLinePointer = fileData.tell()	# Store the pointer to this line as it will be needed later
	breakFlag = True
	while breakFlag == True:
		line = readtabline(fileData)
		if line[0] == '-999':			# Read lines until the -999 string is detected
			endDataCardLineNumber = lineNumber(fileData)	# Store the line number of the end data card
			nFault = endDataCardLineNumber - startDataCardLineNumber - 1	# Compute number of faults
			breakFlag  = False

	faultData = np.zeros((nFault,5))

	fileData.seek(startDataCardLinePointer)	# Go back to the start data card line

	for i in range(nFault):
		line = readtabline(fileData)
		faultData[i,0] = int(line[0])-1
		faultData[i,1] = int(line[1])-1
		faultData[i,2] = int(line[2])-1
		faultData[i,3] = float(line[3])
		faultData[i,4] = float(line[4])

	print(' Done.')

# END OF FILE MANIPULATION: CLOSING FILE --------

# (3.2) ACQUIRING POWER MEASUREMENTS {{{2
if measurementsFile:
	with open(measurementsFile,'r') as fileData:

		breakFlag = True

		# (3.2.1) Searching for start card

		while breakFlag==True:
			line = readtabline(fileData)
			if line[0] == 'BEGGINING ACTIVE POWER MEASURES': breakFlag = False

		print(' --> Beginning active power measurements data reading... ',end='' )

		fileData.readline()
		isP = np.zeros((nBus,nBus),dtype=int)
		wP = np.zeros((nBus,nBus))
		P = np.zeros((nBus,nBus))

		# (3.2.2) Acquiring active power measures
		breakFlag = True
		line = readtabline(fileData)
		while breakFlag == True:
			fromBus = int(float(line[0]))
			toBus = int(float(line[1]))
			isP[fromBus-1,toBus-1] = 1
			P[fromBus-1,toBus-1] = float(line[2])
			wP[fromBus-1,toBus-1] = float(line[3])

			line = readtabline(fileData)
			if line[0] == '-999': breakFlag  = False

		# (3.2.3) Acquiring reactive power measures
		print(' Done.\n --> Beginning reactive power measurements data reading... ',end='' )

		breakFlag = True
		while breakFlag==True:
			line = readtabline(fileData)
			if line[0] == 'BEGGINING REACTIVE POWER MEASURES': breakFlag = False

		wQ = np.zeros((nBus,nBus))
		Q = np.zeros((nBus,nBus))

		fileData.readline()
		breakFlag = True
		line = readtabline(fileData)
		while breakFlag == True:

			fromBus = int(float(line[0]))
			toBus = int(float(line[1]))
			Q[fromBus-1,toBus-1] = float(line[2])
			wQ[fromBus-1,toBus-1] = float(line[3])	

			line = readtabline(fileData)
			if line[0] == '-999':
				breakFlag  = False

		# (3.2.4) Acquiring voltage measures
		print(' Done.\n --> Beginning voltage measurements data reading... ',end='' )

		breakFlag = True
		while breakFlag==True:
			line = readtabline(fileData)
			if line[0] == 'BEGGINING VOLTAGE MEASURES': breakFlag = False
		fileData.readline()

		isV = np.zeros(nBus,dtype=int)
		wV = np.zeros(nBus)
		Vm = np.zeros(nBus)
		breakFlag = True
		line = readtabline(fileData)
		while breakFlag == True:

			fromBus = int(float(line[0]))
			isV[fromBus-1] = 1
			Vm[fromBus-1] = float(line[1])
			wV[fromBus-1] = float(line[2])

			line = readtabline(fileData)
			if line[0] == '-999': breakFlag  = False	

		# (3.2.5) Acquiring available pseudo-measures
		if obs:
			print(' Done.\n --> Beginning available pseudo-measurements data reading... ',end='' )

			breakFlag = True
			while breakFlag==True:
				line = readtabline(fileData)
				if line[0] == 'BEGGINING AVAILABLE PSEUDO-MEASURES (IN ORDER OF PRIORITY)': breakFlag = False
			fileData.readline()

			# Preallocating isPseudo with one row which will be deleted after
			isPseudo = np.zeros((1,2))
			pseudoP = np.zeros((nBus,nBus))
			pseudoWP = np.zeros((nBus,nBus))
			pseudoQ = np.zeros((nBus,nBus))
			pseudoWQ = np.zeros((nBus,nBus))

			breakFlag = True
			while breakFlag == True:
				line = readtabline(fileData)
				if line[0] == '-999': breakFlag  = False

				if breakFlag == True:
					fromBus = int(float(line[0]))
					toBus = int(float(line[1]))
					isPseudo = conc((isPseudo,array([[fromBus-1,toBus-1]])),axis=0)
					pseudoP[fromBus-1,toBus-1] = float(line[2])
					pseudoWP[fromBus-1,toBus-1] = float(line[3])	
					pseudoQ[fromBus-1,toBus-1] = float(line[4])
					pseudoWQ[fromBus-1,toBus-1] = float(line[5])

			# Deleting first row from preallocating
			isPseudo = np.delete(isPseudo,0,axis=0)


		#(3.2.6) Counting P and V readings
		numP = isP.sum() # Number of P readings
		numV = isV.sum()

		print(' Done.\n')

# END OF FILE MANIPULATION: CLOSING FILE --------
# -------------------------------------------------
# (4) TESTING OBSERVABILITY {{{1
# -------------------------------------------------

if obs:
	# (3.3.1) Building Hdelta and gqain matrix
	if verbose > 0: print('\n' + '-'*50 + '\n OBSERVABILITY ANALYSIS \n' + '-'*50)
	else: print('\n --> Beggining observability test.',end='')

	obs = np.zeros((numP,nBus),dtype=float)

	i = 0
	for m in range(nBus):
		for n in range(nBus):
			if isP[m,n]:
				if m == n:
					for k in range(nBus):
						if k != m and k in K[m]: obs[i,k] = -1
						elif k == m: obs[i,k] = len(K[m])
				else:
					for k in range(nBus):
						if k==m: obs[i,k] = 1
						if k==n: obs[i,k] = -1
				i += 1

	gain = transpose(obs)@obs

	# (3.3.2) Printing results
	if verbose > 2: print('\n --> H matrix: H = \n{0}'.format(obs))

	Pt,Lt,Ut = LU(obs)
	dU = array([Ut[k,k] for k in range(Ut.shape[0])])

	if verbose > 2: print('\n --> Upper-triangularized linear-model mF.Jacobian H: HTri = \n{0}'.format(Ut))
	#print('\n --> Main diagonal of U: diag(U) = \n{0}'.format(dU))

	if verbose > 2: print('\n --> Gain matrix: GTri = \n{0}'.format(gain))

	Pt,Lt,Ut = LU(gain)
	dU = array([Ut[k,k] for k in range(Ut.shape[0])])

	if verbose > 1: print('\n --> Upper-triangular factorization of G: U = \n{0}'.format(Ut))
	#print('\n --> Main diagonal of U: diag(U) = \n{0}'.format(dU))

	# (3.3.3) Checking for null pivots and restoring observability from avalable pseudos
	zeroPivots = 0	# Number of null pivots: should be one in order for the network to be observable
	for k in range(dU.shape[0]):
		if abs(dU[k]) < 1e-8:	# Due to numerical error propagation, the null pivots do not come as zero per se, but small numbers. This threshold of 1e-8 was empirically set and seems to work fine for all cases.
			zeroPivots += 1

	if zeroPivots == 1:
		print('\n --> The system seems to be observable: single null pivot detected.')
	else:
		# Restoring observability from available pseudo-measures
		pause('\n --> The system appears not to be observable: there were {0} null pivots. Attempting to restore observability from available pseudo-measures.'.format(zeroPivots))

		# Copy isP
		tempIsP = copy(isP)
		while zeroPivots != 1:
			for pseudoTestCounter in range(isPseudo.shape[0]):
				fromBus = int(isPseudo[pseudoTestCounter,0])
				toBus = int(isPseudo[pseudoTestCounter,1])
				
				# Adding pseudo-measure to the measures roster
				tempIsP[fromBus,toBus] = 1

				try:
					addedPseudoNumber = addedPseudos.shape[1]
				except NameError:
					addedPseudoNumber = 0

				# Building new observability and gain matrices
				obs = np.zeros((numP + addedPseudoNumber + 1,nBus))

				i = 0
				for m in range(nBus):
					for n in range(nBus):
						if tempIsP[m,n]:
							if m == n:
								for k in range(nBus):
									if k != m and k in K[m]: obs[i,k] = -1
									elif k == m: obs[i,k] = len(K[m])
							else:
								for k in range(nBus):
									if k==m: obs[i,k] = 1
									if k==n: obs[i,k] = -1
							i += 1

				gain = transpose(obs)@obs

				Pt,Lt,Ut = LU(gain)
				dU = array([Ut[k,k] for k in range(Ut.shape[0])])

				# Checking for null pivots
				tempZeroPivots = 0
				for k in range(dU.shape[0]):
					if abs(dU[k]) < 1e-8:
						tempZeroPivots += 1

				# Now, if the new pseudo-measure added a new information, add this measure to the set
				if tempZeroPivots < zeroPivots:
					zeroPivots = tempZeroPivots
					if verbose > 0: print(' --> Measure ({0:3},{1:3}) added a new non-null pivot. This measure is added to the measures roster.'.format(fromBus+1,toBus+1))
					try:
						addedPseudos = conc((addedPseudos,array([[fromBus,toBus]])),axis=0)
					except NameError:
						addedPseudos = array([[fromBus,toBus]])
				# But if it did not add anu information, remove this measure
				else: 
					tempIsP[fromBus,toBus] = 0
					if verbose > 0: print(' --> Measure ({0:3},{1:3}) did not add new non-null pivot. Discarding this measure.'.format(fromBus+1,toBus+1))
					if pseudoTestCounter == isPseudo.shape[0]:
						sys.exit('--> Fatal error: the available pseudo-measures were not able to restore observability!')

		# (3.3.4) Updating system measurement parameters and printing added pseudos
		# Updating P, Q and V measurements with pseudos
		isP = copy(tempIsP)

		P += pseudoP
		wP += pseudoWP
		Q += pseudoQ
		wQ += pseudoWQ
		numP = isP.sum() # Number of P readings
		numV = isV.sum()

		print('\n --> Summary: added pseudo-measures: \n{0}'.format(addedPseudos + np.ones((addedPseudos.shape[0],addedPseudos.shape[1]))))

# -------------------------------------------------
# (5) TESTING FOR CRITICAL MEASUREMENTS {{{1
# -------------------------------------------------

if crit:
	if verbose > 0: print('\n' + '-'*50 + '\n TESTING FOR CRITICAL MEASURES \n' + '-'*50)
	else: print('\n --> Testing for critical measures.', end = '')

	tempIsP = copy(isP)

	# Sweeping through available measurements
	for p in range(nBus):
		for q in range(nBus):
			if isP[p,q]:
				if verbose > 0: print(' --> Removing measure ({0:3},{1:3})...'.format(p+1,q+1), end='')
				# Temporarily removing measure (p,q)
				tempIsP[p,q] = 0

				# Building observability matrix again
				obs = np.zeros((numP-1,nBus),dtype=float)

				i = 0
				for m in range(nBus):
					for n in range(nBus):
						if tempIsP[m,n]:
							if m == n:
								for k in range(nBus):
									if k != m and k in K[m]: obs[i,k] = -1
									elif k == m: obs[i,k] = len(K[m])
							else:
								for k in range(nBus):
									if k==m: obs[i,k] = 1
									if k==n: obs[i,k] = -1
							i += 1

				# Calculating gain matrix and triangularizing
				gain = transpose(obs)@obs
				Pt,Lt,Ut = LU(gain)
				dU = array([Ut[k,k] for k in range(Ut.shape[0])])

				# Testing for null pivots
				zeroPivots = 0
				for k in range(dU.shape[0]):
					if abs(dU[k]) < 1e-8:
						zeroPivots += 1

				if zeroPivots > 1:
					print('\n  !-> The system has {2:3} null pivots. Measure ({0:3},{1:3}) is critical.\n'.format(p+1,q+1,zeroPivots))

					# Try to add this measure to the critical measure roster
					try:
						criticalMeas = conc((criticalMeas,array([[p,q]])),axis=0)
					# In case the array criticalMeas does not exist, catch exception and create it
					except NameError:
						criticalMeas = array([[p,q]])

				elif zeroPivots == 1:
					if verbose > 0: 
						print(' Gain matrix has exactly one null pivot.')
				elif zeroPivots == 0: print('\n\n\n ERROR: zero null pivots detected!')

				tempIsP[p,q] = 1

	# Try to print critical measurement list
	try:
		print('\n --> Summary: list of critical measurements: \n{0}'.format(criticalMeas + np.ones((criticalMeas.shape[0],criticalMeas.shape[1]))))
	# In case no critical measurement was found, then this will throw an error
	except NameError:
		print('\n --> No critical measurements were found.')

# -------------------------------------------------
# (7) SETTING UP NUMERICAL METHOD FOR STATE ESTIATION {{{1
# -------------------------------------------------
if estimatestate:

	# (3.5.2) DESCRIBING MEASUREMENT ARRAY Z = [P,Q,V].

	Zflat = mF.Z(P,Q,isP,Vm,isV)

	if verbose > 1: 
		pause('-'*50 + '\n DESCRIBING MEASURING ARRAY\n' + '-'*50)
		for j in range(nBus):	
			for n in range(nBus):
				if isP[j,n]:
					print('Z[{0:2.0f}] = P[{1:3},{2:3}] = {3:>10f}'.format(i+1,j+1,n+1,P[j,n]))
					i+=1
		for j in range(nBus):	
			for n in range(nBus):
				if isP[j,n]:
					if verbose > 1: print('Z[{0:2.0f}] = Q[{1:3},{2:3}] = {3:>10f}'.format(i+1,j+1,n+1,Q[j,n]))
					i += 1

		for j in range(nBus):	
			if isV[j]:
				if verbose > 1: print('Z[{0:2.0f}] = V[{1:3}]     = {2:>10f}'.format(i+1,j+1,Vm[j]))
				i += 1

	# (3.5.3) PRINTING RESULTS
	if verbose > 1:
		print('\n' + '-'*50 + '\n DESCRIBING SYSTEM \n' + '-'*50)
		print(' --> r = {0}\n\n --> x = {1}\n\n --> g = \n{2}\n\n --> b = \n{3}\n\n --> bsh = \n{4}\n\n --> Y = \n{5} \n\n --> stdP = \n{6}\n\n --> stdQ = \n{7}\n\n --> a = \n{8}\n'.format(r,x,g,b,bsh,Y,wP,wQ,a))

	# (3.5.1) CONSTRUCTING W MATRIX
	W = np.empty(2*numP + numV)			# Preallocating
	i = 0						# Starting counter
	for j in range(nBus):	
		for n in range(nBus):
			if isP[j,n]:
				W[i] = 1/wP[j,n]**2
				i+=1
	for j in range(nBus):	
		for n in range(nBus):
			if isP[j,n]:
				W[i] = 1/wQ[j,n]**2
				i += 1

	for j in range(nBus):	
		if isV[j]:
			W[i] = 1/wV[j]**2
			i += 1

	W = np.diag(W)

	# (5.1) Initial guess: flat start
	V = np.ones(nBus)
	theta = np.zeros(nBus)
	H = mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV)

	if verbose > 1:
		print('\n' + '-'*50 + '\n STATE ESTIMATION FLAT START \n' + '-'*50)
	else:
		print(' --> Beggining flat start calculations... ', end='')

	if verbose > 1:
		print('\n --> J = \n{0}\n\n --> h = \n{1}'.format(H,mF.h(V,theta,K,a,y,Y,bsh,isP,isV)))

	# (5.2) Defining simulation parameters from the arguments passed to the script
	absTol = args.tol
	deltaMax = args.deltamax
	maxIter = args.maxiter

	# (5.3) Initiating iteration counter
	itCount = 0	

	# (5.4) Calculating mF.Jacobian and its gradient at flat start

	grad = -transpose(H) @ W @ (mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))

	if verbose <= 1: print('Done.')

	# (5.5) Printing method parameters

	if verbose > 1:
		print('\n --> Norm of the mF.Jacobian gradient at flat start: |grad(J)| = {0}'.format(norm(grad)))

	if verbose > 0:
		pause('\n --> Next step is the power flow numerical method. Press <ENTER> to continue.')
		print('\n' + '-'*50 + '\n BEGGINING STATE ESTIMATION NUMERICAL METHOD \n' + '-'*50)
		print('\n --> Newton-Raphson method parameters: tolerance = {0}, maximum step = {1}, maximum interations = {2}'.format(absTol,deltaMax,maxIter))
	else:
		print('\n --> Beggining numerical method... ', end='')

	# (5.6) Start time counte6
	tstart = time.time()

	# -------------------------------------------------
	# (8) STARTING NUMERICAL METHOD FOR STATE ESTIMATION{{{1
	# -------------------------------------------------

	while(True):
	# (6) STARTING ITERATIONS

		# (6.1) Increasing iteration counter
		itCount += 1
		
		# (6.2) Printing iteration report
		if verbose > 0: print('\n ==> Iteration #{0:3.0f} '.format(itCount) + '-'*50)
		
		# (6.3) Calculating mF.Jacobian
		H = mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV)
		
		# (6.4) Calculating gain matrix
		U = transpose(H) @ W @ H

		# (6.5) Calculating state update
		dX = inv(U) @ transpose(H) @ W @ (mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))

		# (6.6) Updating V and theta
		theta[1:] += dX[0:nBus-1,0]
		V += dX[nBus-1:,0]

		# (6.7) Calculating mF.Jacobian gradient
		grad = -transpose(H) @ W @ (mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))
		
		# (6.8) Printing iteration results
		if verbose > 0:
			print('\n --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}\n --> |grad(J)| = {3}'.format(norm(dX),transpose(dX[nBus-1:,0]),dX[0:nBus-1,0],transpose(norm(grad))))
		if verbose > 1:
			print('\n --> J = \n{0},\n\n --> r = Z - h = \n{2}*\n{1}'.format(mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV),(mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))/norm(mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV)),norm(mF.Z(P,Q,isP,Vm,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))))

		# (6.9) Testing for iteration sucess or divergence
		if norm(dX) < absTol: # (6.8.1) If success
			print('\n --> Sucess!')
			print('\n' + '-'*50 + '\n STATE ESTIMATION CONVERGENCE RESULTS \n' + '-'*50)
			print('\n --> Numerical method stopped at iteration {1} with |deltaX| = {0}.\n --> |grad(J)| = {2}'.format(norm(dX),itCount,norm(grad)))
			break	

		if norm(dX) > deltaMax: #(6.8.2) If diverted
			print('\n --> Error: the solution appears to have diverged on iteration number {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break
	
		if itCount > maxIter - 1: # (6.8.3) If reached maximum iteration
			print('\n --> Error: the solution appears to have taken too long. Final iteration: {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break

	# (6.10) Pausing for each iteration
	if verbose>1:
		pause('\n --> Program paused for next iteration. ')


	# (6.11) Calculating elapsed time
	elapsed = time.time() - tstart

	# (6.12) Calculating normalized residual
	r = mF.Z(P,Q,isP,np.zeros(nBus),isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV)
	cov = inv(W) - H @ inv(U) @ transpose(H)
	rn = array([ r[i]/sqrt(cov[i,i]) for i in range(len(r))])

	# (6.13) Printing convergence results
	print(' --> Final result:\n\n	theta = {0} \n\n	V     = {1}'.format(theta,V))

	if verbose > 1: print('\n --> Residual = \n{2}\n\n --> Normalized residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,rn,r))

#--------------------------------------------------
# (9) SETTING UP NUMERICAL METHOD FOR POWER FLOW {{{1
# -------------------------------------------------
if powerflow or dyn:

	print('\n' + '-'*50 + '\n POWER FLOW METHOD \n' + '-'*50)
	# (5.1) Initial guess: flat start
	V = .9*np.ones(nBus)
	theta = np.zeros(nBus)
	P = np.diag(Pger - Pload)
	Q = np.diag(Qger - Qload)
	isP = np.eye(nBus)
	isV = np.zeros(nBus)

	numP = int(isP.sum())
	numV = int(isV.sum())

	# (3.5.2) DESCRIBING MEASUREMENT ARRAY Z = [P,Q,V].
	if verbose > 1: pause('-'*50 + '\n DESCRIBING MEASURING ARRAY\n' + '-'*50)

	if verbose > 1:
		print('\n' + '-'*50 + '\n POWER FLOW FLAT START \n' + '-'*50)
	else:
		print('\n --> Beggining power flow calculations... ', end='')

	if verbose > 1:
		print('\n --> J = \n{0}\n\n --> h = \n{1}'.format(mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV),mF.h(V,theta,K,a,y,Y,bsh,isP,isV)))

	# (5.2) Defining simulation parameters from the arguments passed to the script
	absTol = args.tol
	deltaMax = args.deltamax
	maxIter = args.maxiter

	# (5.3) Initiating iteration counter
	itCount = 0	

	if verbose > 0:
		pause('\n --> Next step is the numerical method. Press <ENTER> to continue.')
		print('\n' + '-'*50 + '\n BEGGINING POWER FLOW NUMERICAL METHOD \n' + '-'*50)
		print('\n --> Newton-Raphson method parameters: tolerance = {0}, maximum step = {1}, maximum interations = {2}'.format(absTol,deltaMax,maxIter))
	else:
		print('\n --> Beggining numerical method... ', end='')

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
		H = mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV)
		H = np.delete(H,0,axis=0) # Removing slack bar angle derivatives for it is known

		# (6.5) Calculating state update
		deltaSLC = np.delete(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV),0,axis=0)
		dX = inv(H) @ deltaSLC

		# (6.6) Updating V and theta
		theta[1:] += dX[0:nBus-1,0]
		V += dX[nBus-1:,0]
		
		# (6.8) Printing iteration results
		if verbose > 0:
			print(' --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}'.format(norm(dX),transpose(dX[nBus-1:,0]),dX[0:nBus-1,0]))
		if verbose > 1:
			print('\n --> J = \n{0},\n\n --> r = Z - h = \n{2}*\n{1}'.format(	mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV),
												(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))/norm(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV)),
												norm(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))))

		# (6.9) Testing for iteration sucess or divergence
		if norm(dX) < absTol: # (6.8.1) If success
			print('\n --> Sucess!')
			print('\n --> Numerical method stopped at iteration {1} with |deltaX| = {0}.'.format(norm(dX),itCount))
			break	

		if norm(dX) > deltaMax: #(6.8.2) If diverted
			print('\n --> Error: the solution appears to have diverged on iteration number {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break
		
		if itCount > maxIter - 1: # (6.8.3) If reached maximum iteration
			print(' --> Error: the solution appears to have taken too long. Final iteration: {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break

		# (6.10) Pausing for each iteration
		if verbose>1:
			pause('\n --> Program paused for next iteration. ')


	# (6.11) Calculating elapsed time
	elapsed = time.time() - tstart

	# (6.13) Printing convergence results
	print(' --> Final result:\n\n	theta = {0} \n\n	V     = {1}'.format(theta,V))

	if verbose > 1: print('\n --> Residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,r))

# -------------------------------------------------
# (11) STARTING NUMERICAL METHOD FOR DYNAMICAL SIMULATION {{{1
# -------------------------------------------------
if dyn:
	# Gathering simulation parameters from args
	tFinal = args.time
	pps = args.pointspersecond

	print('\n' + '-'*50 + '\n DYNAMIC FAULT SIMULATION \n' + '-'*50)
	print('\n --> Setting up dynamical simulation...')

	# (11.1) Pre-fault state
	YLoad = [(Pload[i] - 1j*Qload[i])/V[i]**2 for i in range(nBus)]	# YLoad is the equivalent conductance load matrix. Each load is modelled as a constant impedance
	
	# Permutating the Y and K matrixes. Explanation: the Y matrix built before did not consider which bars were generators and which are not. Nevertheless, to build the Y_RED matrix, the generator buses must be the first rows. This means that we must obtain a new Y matrix wherein the first n columns are the genertor buses.
	# In order to do this, a permutation matrix PDyn is built. The building process works as follows:
	PDyn = np.eye(nBus)
		
	genDataDyn = copy(genData)
	busNameDyn = copy(busName)
	for i in range(nGen):
		for k in range(int(genDataDyn[i,0])):
			if not mF.isGen(k,genDataDyn):
				temp = copy(PDyn[int(genDataDyn[i,0]-1)])
				PDyn[int(genDataDyn[i,0]-1)] = copy(PDyn[k])
				PDyn[k] = temp
				genDataDyn[i,0] = k+1
				break

	# Remember that we need to permute both columns and rows. To permute rows, it is needed to left-multiply Y by P (PY) and to columns swap, right-multiply (YP).
	YDyn = PDyn@Y@PDyn
	YLoadDyn = PDyn@YLoad
	VDyn,thetaDyn = PDyn@V, PDyn@theta

	Yred,C,D = mF.reduceGrid(YDyn,YLoadDyn,VDyn,genData)

	
	# Pre-fault mechanical power
	pm = np.zeros(nGen)
	for i in range(nGen):
		pm[i] = (VDyn[i]**2)*real(Yred[i,i]) + sum( [ (C[i,j]*sin(thetaDyn[i] - thetaDyn[j]) + D[i,j]*cos(thetaDyn[i] - thetaDyn[j])) for j in range(nGen)])

	fig = [ [] for i in range(nFault)]

	for q in range(nFault):

		faultData[q,4] /= 1000

		lineAdmittance = Y[int(faultData[q,1]),int(faultData[q,2])] 

		YFault = copy(Y)
		YLoadFault = copy(YLoad)
	
		YFault[int(faultData[q,1]),int(faultData[q,2])] = 0
		YFault[int(faultData[q,2]),int(faultData[q,1])] = 0

		if int(faultData[q,3]) == 0:
			YLoadFault[int(faultData[q,1])] += lineAdmittance*1e3
			YLoadFault[int(faultData[q,2])] += lineAdmittance
		elif int(faultData[q,3]) == 1:
			YLoadFault[int(faultData[q,1])] += lineAdmittance
			YLoadFault[int(faultData[q,2])] += lineAdmittance*1e3
		else:
			YLoadFault[int(faultData[q,1])] += lineAdmittance/faultData[q,3]
			YLoadFault[int(faultData[q,2])] += lineAdmittance/(1-faultData[q,3])
		YDynFault = PDyn@YFault@PDyn
		YLoadDynFault = PDyn@YLoadFault
		VDynFault,thetaDynFault = PDyn@V, PDyn@theta

		YredFault,CFault,DFault = mF.reduceGrid(YDynFault,YLoadDynFault,VDynFault,genDataDyn)

		print(' --> Starting numerical simulation of fault {0}...'.format(int(faultData[q,0] + 1)),end='')


		tFault = np.linspace(0,faultData[q,4],int(round(pps*(faultData[q,4]))),endpoint=True)
		y0 = np.zeros(2*nGen)

		for i in range(nGen): y0[2*i] = thetaDyn[i]

		solFault = odeint(mF.sm2,y0,tFault,args=(CFault,DFault,YredFault,VDyn,pm,genDataDyn))

		# Calculating post-fault system

		YPost = copy(YDyn)
		YLoadPost = copy(YLoadDyn)

		YPost[int(faultData[q,1]),int(faultData[q,2])] = 0
		YPost[int(faultData[q,2]),int(faultData[q,1])] = 0

		YDynPost = PDyn@YFault@PDyn
		YLoadDynPost = PDyn@YLoadPost

		YredPost,CPost,DPost = mF.reduceGrid(YDynPost,YLoadDynPost,VDyn,genDataDyn)

		tPost = np.linspace(faultData[q,4],tFinal,int(round(pps*(tFinal - faultData[q,4]))),endpoint=True)
		solPost = odeint(mF.sm2,solFault[-1,:],tPost,args=(CPost,DPost,YredPost,VDyn,pm,genDataDyn))

		print('Done.')

		# Plotting results

		sol = conc((solFault,solPost),axis=0)
		t = conc((tFault,tPost),axis=0)

		fig[q] = pyplt.figure()
		ax1 = fig[q].add_subplot(2,1,1)
		ax2 = fig[q].add_subplot(2,1,2)

		for i in range(nGen):
			ax2.plot(t,sol[:, 2*i])
			ax1.plot(t,sol[:, 2*i+1]-sol[:, 1],label=busName[int(genData[i,0])-1])

		ax1.axvline(x=faultData[q,4], color='k', lw=1, linestyle='--',label=r'$t_F$'.format(faultData[q,4]))
		ax2.axvline(x=faultData[q,4], color='k', lw=1, linestyle='--',label=r'$t_F$'.format(faultData[q,4]))
		ax1.legend(loc="lower right")

		ax1.set_ylabel(r'Synchronized angular speed $\left(\sfrac{rad}{s}\right)$')
		ax2.set_ylabel(r'Power angle $\left(rad\right)$')
		ax1.grid(which='major',axis='both')
		ax1.grid(which='major',axis='both')

		ax1.grid(which='major',axis='both')
		ax2.grid(which='major',axis='both')

		fig[q].suptitle(r'Fault {0} results. Description: short-circuit at line connecting buses {1} and {2}, at position {3}, with open time {4} ms'.format(int(faultData[q,0])+1, int(faultData[q,1]+1), int(faultData[q,2]+1), faultData[q,3], 1000*faultData[q,4]))
		ax2.set_xlabel('Time (s)')

	print(' --> Plotting dynamical simulation results.')
	pyplt.legend()

	pyplt.show()
