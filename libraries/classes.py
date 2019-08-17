# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: classes library for pyEPSim
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.0
# DESCRIPTION: this file describes the classes used in the pyEPSim program.
# -------------------------------------------------

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

from copy import deepcopy
copy = deepcopy

from libraries.tabulate.tabulate import tabulate

# (1) Case object {{{1
# The case object stores data about a particular "case". "Cases" are a strain of data which comprise simulation parameters of a given system. These parameters include topological data like bus, branch and generator data and simulatory parameters, such as simulation time, plot data, tolerances et cetera. This object serves the purpose of giving the user an easy way to modify and customize the set of parameters -- be them topological or simulatory -- without having to generate a whole new file. In general, a case file is built from a *.net supplied by the user; see the loadCase function in file loadCase.py
# --> "name" is a human-readable name that can be remembered by the user;
# --> "busData" is a list of "bus" objects (see (2)), which store information about buses like injected power, load power, name, number;
# --> "branchData" is a list of "branch" objects (see (3)) which store information about the branches of the system, such as resistance and reactance, transformer turns ratio and so on;
# --> "genData" is a list of "generator" objects (see (4)) which sotre information about system generators. These generators are modelled as synchronous machines; the "generator" object sotres basically constructive parameters of the machines.
# --> "faultData" is a list of "fault" objects (see (5)) which store information about possible faults in the system. The listed faults will be used for dynamical simulation.
# These quantities are defined in the class creation; that is, when a case instance is created, these parameters need to be informed. In general, these quantities will be declared in a net file (*.net) and the file is loaded through the loadCase funtion in loadCase.py . Apart from these values, there are other constructive/topological parameters of the system that must be taken in place. The numerical method uses matrixes and not classes to perform calculations; these matrixes have to be constructed from the case parameters. The matrixe are explained below.
# --> z is a matrix defining the branch impedances; z[m,n] is the branch impedance from bus m to bus n; if buses m or n are not connected, then z[m,n] is zero. Also, z has a null main diagonal.
# --> y is the element-wise inverse of z; g being its real part, b its imaginary.
# --> bsh is the matrix of shunt capacitances between buses; bsh[m,n] is the shunt capacitance of line between buses m and n (zero if m and n are not neighbors); bsh[k,k] is the shunt capacitance of bus k.
# --> Y is the multi-dimensional thevenin equivalent susceptance of the grid.

class case:
	def __init__(self,name, busData,branchData,genData,faultData,Sb,Vb):
		self.name = name
		self.busData = busData
		self.branchData = branchData
		self.genData = genData
		self.faultData = faultData	
		self.Sb = Sb			# Base power value
		self.Vb = Vb			# Base voltage value

		self.nBus = len(self.busData)

		self.r, self.x, self.a, self.bsh, self.K  = self.updateMatrixes()		
		self.z = array([[self.r[m,n] + 1j*self.x[m,n] for m in range(self.nBus)] for n in range(self.nBus)])
		self.y = array([[1/self.z[m,n] if self.z[m,n] != 0 else 0 for m in range(self.nBus)] for n in range(self.nBus)])
		self.g = real(copy(self.y))
		self.b = imag(copy(self.y))
		self.bsh = array([[ self.bsh[m,n]/2 if m!=n else self.bsh[m,n] for n in range(self.nBus)] for m in range(self.nBus)])

		self.Y = self.buildY();
		self.B = imag(copy(self.Y))
		self.G = real(copy(self.Y))

	def __str__(self,**kwargs):
		if ('verbose' in kwargs): verbose = kwargs['absTol']
		else: verbose = 0

		# tableFormat is the string that dictates the final table output format. It is derived from the formats of the dependency python-tabulate. See the documentation at the "table format" chapter.		
		if ('tableformat' in kwargs): tFormat = kwargs['tableformat']
		else: tableformat = 'psql' 

		print(' --> Printing case \'{0}\' data:'.format(self.name))
		
		# PRINTING BUS DATA ------------------------
		print('\n >> Bus list')
		tabRows = []
		tabHeader = ['Number', 'Name', 'Active Load (pLoad)', 'Reactive Load (qLoad)', 'Active Generation (pGer)', 'Reactive Generation (qGer)', 'Shunt capacitance (bsh)', 'Shunt conductance (gsh)']
		for bus in self.busData:
			tabRows.append([bus.number, bus.name, bus.pLoad, bus.qLoad, bus.pGen, bus.qGen, bus.bsh, bus.gsh ])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

		# PRINTING BRANCH DATA ---------------------
		print('\n\n >> Branch list')
		tabRows = []
		tabHeader = ['From Bus', 'To Bus', 'Resistance (r)', 'Reactance (x)', 'Shunt susceptance (bsh)', 'Transformer turns ratio (a)']
		for branch in self.branchData:
			tabRows.append([branch.fromBus, branch.toBus, branch.r, branch.x, branch.bsh, branch.a])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

		# PRINTING GENERATOR DATA ------------------
		print('\n\n >> Generator list')
		tabRows = []
		tabHeader = ['Bus Number', 'Rated Power', 'H', 'D', 'ra', 'xL', 'xd', 'xPd', 'xPPd', 'tPdo', 'tPPdo', 'xq', 'xPq', 'xPPq', 'tPqo', 'tPPqo']
		for gen in self.genData:
			tabRows.append([gen.busNumber, gen.ratedPower,gen.H, gen.D, gen.ra, gen.xL, gen.xd, gen.xPd, gen.xPPd, gen.tPdo, gen.tPPdo, gen.xq, gen.xPq, gen.xPPq, gen.tPqo, gen.tPPqo])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))	

		# PRINTING FAULT DATA ----------------------
		print('\n\n >> Possible faults list')
		tabRows = []
		tabHeader = ['Branch Number', 'Location', 'Opening Time']
		for fault in self.faultData:
			tabRows.append([fault.branch, fault.location, fault.openingTime])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

		return ''


# Function case.updateMatrixes is meant to update matrixes r,x,a,bsh and K everytime these are called, or everytime they are needed. This is done to prevent inconsistencies if the user changes a variable directly, say for example:
# case.branchData[1].r = 5
# In this case the user has directly updated the value of a branch resistance; all matrixes must be re-calculated. Function updateMatrixes is called whenever case.r, case.x, case.a, case.bsh, case.K are called, so they are recalculated for the present values of bus and branch data.
	def updateMatrixes(self):
		r = np.zeros((self.nBus,self.nBus))
		x = np.zeros((self.nBus,self.nBus))
		a = np.ones((self.nBus,self.nBus))
		bsh = np.zeros((self.nBus,self.nBus))
		K = [ [] for i in range(self.nBus)]
		
		for branch in self.branchData:
			K[branch.fromBus-1] += [branch.toBus - 1]
			K[branch.toBus-1] += [branch.fromBus - 1]

			r[branch.fromBus-1,branch.toBus-1] = branch.r
			r[branch.toBus-1,branch.fromBus-1] = branch.r

			x[branch.fromBus-1,branch.toBus-1] = branch.x
			x[branch.toBus-1,branch.fromBus-1] = branch.x

			bsh[branch.fromBus-1,branch.toBus-1] = branch.bsh	
			bsh[branch.toBus-1,branch.fromBus-1] = branch.bsh

			if branch.a != 0: a[branch.fromBus-1,branch.toBus-1] = branch.a

		for k in range(self.nBus):
			bsh[k,k] = self.busData[k].bsh

		return r,x,a,bsh,K

	def buildY(self):
		Y = -self.a*transpose(self.a)*self.y
		for k in range(self.nBus): Y[k,k] = 1j*self.bsh[k,k] + sum([self.a[k,m]**2*self.y[k,m] + 1j*self.bsh[k,m] for m in self.K[k]])
		return Y		

# (2) Bus object {{{1
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

# (3) Branch object {{{1
# The branch object stores data for a particular branch of the system:
#--> "fromBus" and "toBus" are the numbers of respectively the first and second buses attached to the branch.
#--> "r"  and "x" are respectively the equivalent resistance and reactance of the branch.
#--> "bsh" is the shunt susceptance of the branch according to the pi model
#--> "a" is the turns ratio of the transformer attached to the branch when there is one.
class branch:
	def __init__(self,fromBus,toBus,r,x,bsh,a):
		self.fromBus = int(fromBus)
		self.toBus = int(toBus)
		self.r = float(r)
		self.x = float(x)
		self.bsh = float(bsh)
		self.a = float(a)

# (4) Generator object {{{1
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
	def __init__(self,busNumber,ratedPower,H,D,ra,xL,xd,xPd,xPPd,tPdo,tPPdo,xq,xPq,xPPq,tPqo,tPPqo):
		self.busNumber = busNumber
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

# (5) Fault object {{{1
class fault:
	def __init__(self,branch,location,openingTime):
		self.branch = branch
		self.location = location
		self.openingTime = openingTime
