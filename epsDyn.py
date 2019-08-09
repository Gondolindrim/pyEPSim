#!/usr/bin/python
# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: pyEPSim, a Python dynamical simulator for Electric Power Systems.
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 5.1
# DESCRIPTION: this program is a state estimator for an electric power system.
# The program is designed to dynamically simulate an Electric Power System provided the descriptive files of the system through integration of its Algebraic-Differential equations.
# BIBLIOGRAPHY:
# [1] 
# [2]
# [3]
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
# (2) CUSTOM FUNCTIONS, ARGUMENTS AND OBJECTS{{{1
# -------------------------------------------------
# (2.1) Set printout options: this allows the program printout to span the whole terminal {{{2
np.set_printoptions(suppress=False,linewidth=999,threshold=999)

# (2.2) Clear  function: clears the terminal. Similar to MATLAB's cls and bash's cls {{{2
def clear(): os.system('cls' if os.name == 'nt' else 'clear')

# (2.3) "pause" function: halts program execution and asks for enter {{{2
def pause(string):
    programPause = input(string)

# (2.4) readtabline function {{{2
# readtabline reads a line in a given file f then arranges it in an array which vectors are the strings
#	in the line that are separated by the tab character '\t'. It also removes the new line '\n' character
#	from the last slot in the line.
def readtabline(f): return f.readline().strip().split('\t')

# (2.5) Arguments and options {{{2
#	Creates arguments:
#	--net [NETWORK FILE] is the target system network file. This string is stored in a 'networkFile' variable. The network file describes the system topologically and parametrically: how many buses, which buses are connected to which, how many generators and how they are connected; the parameters of these generators, the faults applicable and their quality.
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
networkFile = args.net

# (2.6) Case object {{{2
# The case object stores data about a particular "case". "Cases" are a strain of data which comprise simulation parameters of a given system. These parameters include topological data like bus, branch and generator data and simulatory parameters, such as simulation time, plot data, tolerances et cetera. This object serves the purpose of giving the user an easy way to modify and customize the set of parameters -- be them topological or simulatory -- without having to generate a whole new file.
# --> "busData" is a list of "bus" objects (see (2.7)), which store information about buses like injected power, load power, name, number;
# --> "branchData" is a list of "branch" objects (see (2.8)) which store information about the branches of the system, such as resistance and reactance, transformer turns ratio and so on;
# --> "genData" is a list of "generator" objects (see (2.9)) which sotre information about system generators. These generators are modelled as synchronous machines; the "generator" object sotres basically constructive parameters of the machines.
# --> "faultData" is a list of "fault" objects (see (2.10)) which store information about possible faults in the system. The listed faults will be used for dynamical simulation.
class case:
	def __init__(self,busData,branchData,genData,faultData,Sb,Vb):
		self.busData = busData
		self.branchData = branchData
		self.genData = genData
		self.faultData = faultData	
		self.Sb = Sb			# Base power value
		self.Vb = Vb			# Base voltage value

# (2.7) Bus object {{{2
# The bus object stores data for a particular bus of the net:
# --> "number" is the bus number in a list. This number can be user-attributed in the net-file,
# but such numbers must be consecutive and have no gaps.
# --> "name" is a human-readable name for that particular bus, so that the user can
# distinguish buses by a name rather than their numbers.
# --> "pLoad" and "qLoad" are respectively the active and reactive power that the bus exports as load.
# These loads are further modelled in the algorithm by a constant impedance.
# --> "pGen" and "qGen" are respectively the active and reactive power injected into the bus.
# --> "gsh" and "bsh" are respectively the shunt conductance and susceptance attached to the bus. These
# will be added to the pLoad and qLoad after these are converted to shunt impedances.
class bus:
	def __init__(self,number,name,pLoad,qLoad,pGen,qGen,gsh,bsh):
		self.number = int(number)
		self.name = str(name)
		self.pLoad = float(pLoad)
		self.qLoad = float(qLoad)
		self.pGen = float(pGen)
		self.qGen = float(qGen)
		self.bsh = float(bsh)
		self.gsh = float(gsh)

# (2.8) Branch object {{{2
# The branch object stores data for a particular branch of the system:
#--> "fromBus" and "toBus" are the numbers of respectively the first and second buses attached to the branch.
#--> "r"  and "x" are respectively the equivalent resistance and reactance of the branch.
#--> "bsh" is the shunt susceptance of the branch according to the pi model
#--> "a" is the turns ratio of the transformer attached to the branch when there is one.
class branch:
	def __init__(self,fromBus,toBus,neighbors,r,x,bsh,a):
		self.fromBus = int(fromBus)
		self.toBus = int(toBus)
		self.r = float(r)
		self.x = float(x)
		self.bsh = float(bsh)
		self.a = float(a)

# (2.6) Generator object {{{2
# The generator object stores the parameters of a given generator. Generators are modelled as synchronous machines:
# --> "H" is the inertia constant of the machine, in p.u.;
# --> "D" is the damping coefficient;
# --> "ra" is the equivalent armature resistance;
# --> "xL" is the rotor magnetizing reactance (aka static reactance)
# --> "xq" and "xd" are the quadrature- and direct-axis static reactances of the rotor;
# --> "xPq" and "xPd" are the quadrature- and direct-axis transient reactances of the rotor;
# --> "xPPq" and "xPPd" are quadrature- and direct-axis sub-transient reactances of the rotor;
# --> "tPqo" and "tPdo" are the rotor quadrature- and direct-axis transient time constants;
# --> "tPPqo" and "tPPdo" are rotor quadrature- and direct-axis sub-transient time constants;
class generator:
	def __init__(self,branchNumber,ratedPower,H,D,ra,xL,xd,xPd,xPPd,tPdo,tPPdo,xq,xPq,xPPq,tPqo,tPPqo):
		self.ratedPower = ratedPower
		self.H = H
		self.D = D
		self.ra = ra
		self.xL = xL
		self.xd = xd
		self.xPd = xPd
		self.xPPd = xPPd
		self.tPdo = tPdo
		self.tPPdo = tPPdo
		self.xq = xq
		self.xPq = xPq
		self.xPPq = xPPq
		self.tPqo = tPqo
		self.tPPqo = tPPqo

# (2.7) Fault object {{{2
class fault:
	def __init__(self,branch,location,openingTime):
		self.branch = branch
		self.location = location
		self.openingTime = openingTime

# -------------------------------------------------
# (3) DESCRIBRING SYSTEM: READING FILE {{{1
# -------------------------------------------------

# loadcase is a function that takes the pointer to a file and returns a case structure (see (2.6)).
# The structure is composed of bus, branch, generator and fault data, comprised of the lists of
# buses, branches, generators and faults. Each of these have a particular set of attributes that
# are explained in their respective object definitions.
def loadCase(fileName):
	printf(' --> Loading case file \'{0}\'...'.format(fileName))
	with open(fileName,'r') as fileData:

		line = []
		while line[0] ~= 'MVA base:': fileData.readline().strip().split('\t')
		Sb = float(line[1])	# Extracting power base value

		while line[0] ~= 'V base:': fileData.readline().strip().split('\t')
		Vb = float(line[1])	# Extracting voltage base value

		# Bus data parameters ------------------------------------------
		# Searching for bus data start card
		while line[0] ~= 'BUS DATA FOLLOWS': line = fileData.readline().strip().split('\t')
		fileData.readline()	# Skip the '---' line		

		# Acquiring bus data
		busList = []
		while line[0] ~= '-999':
			line = readtabline(fileData)
			busList.append(bus(line[1],line[7]/Sb,line[8]/Sb,line[9]/Sb,line[10]/Sb,line[15],line[16]))

		# Branch data parameters ---------------------------------------
		# Searching for branch data start card
		while line[0] ~= 'BRANCH DATA FOLLOWS':	line = fileData.readline().strip().split('\t')
		fileData.readline()	# Skip the '---' line	

		# Acquiring branch data
		branchList = []	
		while line[0] ~= '-999': branchList.append(branch(line[0],line[1],line[6],line[7],line[8],1/float(line[14])))

		# Generator data parameters ------------------------------------
		# Searching for generator data start card
		while line[0]~= 'GENERATOR DATA FOLLOWS': line = fileData.readline().strip().split('\t')
		fileData.readline()	# Skip the '---' line	

		genList = []	# genList was named like so because there already is a genData variable in the program
		while line[0] ~= '-999':
			line = fileData.readline().strip().split('\t')
			genList.append(generator(line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]))
		
		# Fault data ---------------------------------------------------
		# Searching for fault data start card 
		while line[0] ~= 'FAULT DATA FOLLOWS': line = readtabline(fileData)
		next(fileData)	# Skip the '---' line	

		faultList = []	# As the same case with genList, faultList was named like so because there already is a faultData variable in the program
		while line[0] ~= '-999':
			line = readtabline(fileData)
			faultList.append(fault( line[0], line[1], line[2])

		print(' Done.')

	return case(busList,branchList,genList,faultList,Sb,Vb)

print('\n' + '-'*50 + '\n READING FILES \n' + '-'*50)

# (3.1.5) Building Y matrix
print('\n --> Building Y matrix... ',end='')

z = array([[r[m,n] + 1j*x[m,n] for m in range(nBus)] for n in range(nBus)])
y = array([[1/z[m,n] if z[m,n] != 0 else 0 for m in range(nBus)] for n in range(nBus)])

g = real(copy(y))
b = imag(copy(y))

bsh = array([[ bsh[m,n]/2 if m!=n else bsh[m,n] for n in range(nBus)] for m in range(nBus)])

Y  = -a*transpose(a)*y

for k in range(nBus): Y[k,k] = 1j*bsh[k,k] + sum([a[k,m]**2*y[k,m] + 1j*bsh[k,m] for m in K[k]])

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

		lineAdmittance = Y[int(faultData[q,1]),int(faultData[q,2])] 

		YFault = copy(Y)
		YLoadFault = copy(YLoad)
		VDynFault = copy(V)
	
		YFault[int(faultData[q,1]),int(faultData[q,2])] = 0
		YFault[int(faultData[q,2]),int(faultData[q,1])] = 0


		if int(faultData[q,3]) == 0:
			#YLoadFault[int(faultData[q,1])] += lineAdmittance*1e3
			YLoadFault[int(faultData[q,1])] = 0
			VDynFault[int(faultData[q,1])] = 0
			YLoadFault[int(faultData[q,2])] += lineAdmittance
		elif int(faultData[q,3]) == 1:
			YLoadFault[int(faultData[q,1])] += lineAdmittance
			YLoadFault[int(faultData[q,2])] += lineAdmittance*1e3
		else:
			YLoadFault[int(faultData[q,1])] += lineAdmittance/faultData[q,3]
			YLoadFault[int(faultData[q,2])] += lineAdmittance/(1-faultData[q,3])

		YDynFault = PDyn@YFault@PDyn
		YLoadDynFault = PDyn@YLoadFault
		VDynFault,thetaDynFault = PDyn@VDynFault, PDyn@theta

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
			ax2.plot(t,sol[:, 2*i] - sol[:, 0])
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
