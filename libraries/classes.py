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

import libraries.powerFlow as pF

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
		self.nGen = len(self.genData)

		self.r, self.x, self.a, self.bsh, self.K  = self.updateMatrixes()		
		self.z = array([[self.r[m,n] + 1j*self.x[m,n] for m in range(self.nBus)] for n in range(self.nBus)])
		self.y = array([[1/self.z[m,n] if self.z[m,n] != 0 else 0 for m in range(self.nBus)] for n in range(self.nBus)])
		self.g = real(copy(self.y))
		self.b = imag(copy(self.y))
		self.bsh = array([[ self.bsh[m,n]/2 if m!=n else self.bsh[m,n] for n in range(self.nBus)] for m in range(self.nBus)])

		self.Y = self.buildY();
		self.B = imag(copy(self.Y))
		self.G = real(copy(self.Y))

	def __str__(self):

		tableformat = 'psql' 

		print(' --> Printing case \'{0}\' data:'.format(self.name))
		
		self.printBusData()
		self.printBranchData()
		self.printGenData()
		self.printFaultData()


		return ''

	# runPowerFlow() runs the power flow of the case and updates the finalVoltage and finalAngle attributes of the bus instances in the busData list. If the power flow method is not successful it returns a false value.
	def runPowerFlow(self):
		V, theta, r, elapsed, itCount, success = pF.powerFlow(self)
		if not success:
			print(' >>> Power flow update for case class instance returned not successful.')
			V, theta = array(['None']*self.nBus), array(['None']*self.nBus),

		for k in range(self.nBus): self.busData[k].finalVoltage, self.busData[k].finalAngle = V[k], theta[k]
		
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
			K[branch.fromBus] += [branch.toBus]
			K[branch.toBus] += [branch.fromBus]

			r[branch.fromBus,branch.toBus] = branch.r
			r[branch.toBus,branch.fromBus] = branch.r

			x[branch.fromBus,branch.toBus] = branch.x
			x[branch.toBus,branch.fromBus] = branch.x

			bsh[branch.fromBus,branch.toBus] = branch.bsh	
			bsh[branch.toBus,branch.fromBus] = branch.bsh

			if branch.a != 0: a[branch.fromBus,branch.toBus] = branch.a

		for k in range(self.nBus):
			bsh[k,k] = self.busData[k].bsh

		return r,x,a,bsh,K

# Method case.buildY builds the admittance matrix Y of the system
	def buildY(self):
		Y = -self.a*transpose(self.a)*self.y
		for k in range(self.nBus): Y[k,k] = 1j*self.bsh[k,k] + sum([self.a[k,m]**2*self.y[k,m] + 1j*self.bsh[k,m] for m in self.K[k]])
		return Y	

# isGen(busN) returns True if the bus with number busN has a generator attached to it.
	def isGen(self,busN):
		for gen in self.genData:
			if int(gen.busNumber) == busN:
				return True

		return False

# The 'reduceMatrixes' method returns the matrixes of the reduced system. In order to do this, the system will need to be reorganized so that the buses attached to a generator are numbered first. Since this method alters the sequence of buses (that is, it reorganizes the number of the buses), it is recommended that this method be run through a copied instance, that is, run 'reducedCase = copy(case)' and then 'reducedCase.reduceGrid()'. Although this may seem like a problem, given that the user will not know which buses were re-numbered, the bus names defined won'be altered, meaning that the human-readable names do not change. This ultimately means that it is imperative to give each bus a human-readable name. 
	def reduceMatrixes(self):#(Y,Yload,V,genData):
		nBus = self.nBus
		nGen = self.nGen

		# Obtaining the bus number permutation matrix through an equivalence matrix. Equiv is a nGen x 2 matrix where n is the generator number. The first column are the 'old' bus numbers and the second is the 'new' bus numbers; that way, each line shows an equivalence between the old and the new bus numbers.
		equiv = array([ i for i in range(nBus)])

		for i in range(nBus):
			for k in range(i):
				if not self.isGen(equiv[k]) and self.isGen(equiv[i]):
					temp = equiv[k]
					equiv[k] = equiv[i]
					equiv[i] = temp
					break

		
		# Building component matrixes

		Y1 = self.Y[0 : nGen, 0 : nGen]
		Y2 = self.Y[0 : nGen , nGen : nBus+1]
		Y3 = self.Y[nGen : nBus , 0 : nGen]
		Y4 = self.Y[nGen : nBus+1 , nGen : nBus+1]

		YLoad = [(bus.pLoad - 1j*bus.qLoad)/bus.finalVoltage**2 for bus in self.busData]	# YLoad is the equivalent conductance load matrix. Each load is modelled as a constant impedance

		Ylg = np.diag(Yload[0:nGen])

		Yll = np.diag(Yload[nGen:nBus])

		Ytrans = np.array([1/(gen.ra + 1j*gen.xPq) for gen in self.genData]) # Y' no livro
		Ytrans = np.diag(Ytrans)
		
		YA = Ytrans

		YB = conc((-Ytrans,np.zeros((nGen,nBus-nGen))),axis=1)
		
		YC = conc((-Ytrans,np.zeros((nBus-nGen,nGen))),axis=0)

		YDtop = conc((Ytrans+Y1+Ylg,Y2),axis=1)
		YDbot = conc((Y3,Y4 + Yll),axis=1)

		YD = conc((YDtop,YDbot),axis=0)

		# Calculating YRED and C and D coefficients

		Yred = YA - YB @ inv(YD) @ YC
		V = array([ bus.finalVoltage for bus in self.busData])
		C = np.array([ [ V[i]*V[j]*imag(Yred[i,j]) if j != i else 0 for j in range(nGen)] for i in range(nGen)])
		D = np.array([ [ V[i]*V[j]*real(Yred[i,j]) if j != i else 0 for j in range(nGen)] for i in range(nGen)])

		return [Yred,C,D]

	def printBusData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' bus list'.format(self.name))
		tabRows = []
		tabHeader = ['Number', 'Name', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)', 'Active Generation\npGer (MW)', 'Reactive Generation\nqGen (MVAR)', 'Shunt capacitance\nbsh (p.u.)', 'Shunt conductance\ngsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (deg)']
		for bus in self.busData:
			tabRows.append([bus.number, bus.name, bus.PVtype, bus.pLoad, bus.qLoad, bus.pGen, bus.qGen, bus.bsh, bus.gsh, bus.finalVoltage, bus.finalAngle*180/np.pi ])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

	def printBranchData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' branch list'.format(self.name))
		tabRows = []
		tabHeader = ['Number','From Bus\n(Tap bus)', 'To Bus\n(Z bus)', 'Resistance\n r (p.u.)', 'Reactance\n x (p.u.)', 'Shunt susceptance\n bsh (p.u.)', 'Transformer\nturns ratio (a)']
		for branch in self.branchData:
			tabRows.append([branch.number,self.busData[branch.fromBus].name, self.busData[branch.toBus].name, branch.r, branch.x, branch.bsh, branch.a])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

	def printGenData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' generator list (sorted by bus number)'.format(self.name))
		tabRows = []
		tabHeader = ['Bus', 'Rated Power\n(MW)', 'H\n(p.u.)', 'D\n(p.u.)', 'ra\n(p.u.)', 'xL\n(p.u.)', 'xd\n(p.u.)', 'xPd\n(p.u.)', 'xPPd\n(p.u.)', 'tPdo\n(s)', 'tPPdo\n(s)', 'xq\n(p.u.)', 'xPq\n(p.u.)', 'xPPq\n(p.u.)', 'tPqo\n(s)', 'tPPqo\n(s)']
		for gen in self.genData:
			tabRows.append([self.busData[gen.busNumber].name, gen.ratedPower,gen.H, gen.D, gen.ra, gen.xL, gen.xd, gen.xPd, gen.xPPd, gen.tPdo, gen.tPPdo, gen.xq, gen.xPq, gen.xPPq, gen.tPqo, gen.tPPqo])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))	

	def printFaultData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql'
 
		print('\n >> Case \'{0}\'  possible faults list'.format(self.name))
		tabRows = []
		tabHeader = ['Branch Number', 'Location', 'Opening Time']
		for fault in self.faultData:
			tabRows.append([fault.branch, fault.location, fault.openingTime])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

# (2) Bus object {{{1
# The bus object stores data for a particular bus of the net:
# --> "number" is the bus number in a list. This number can be user-attributed in the net-file,
# but such numbers must be consecutive and have no gaps.
# --> "name" is a human-readable name for that particular bus, so that the user can
# distinguish buses by a name rather than their numbers.
# --> "type" is thetype of the bus, meaning PV (hold voltage in voltage limits), PQ (hold generation within limits), VT (hold phase), UN (unregulated)
# --> "pLoad" and "qLoad" are respectively the active and reactive power that the bus exports as load.
# These loads are further modelled in the algorithm by a constant impedance.
# --> "pGen" and "qGen" are respectively the active and reactive power injected into the bus.
# --> "gsh" and "bsh" are respectively the shunt conductance and susceptance attached to the bus. These
# will be added to the pLoad and qLoad after these are converted to shunt impedances.
# --> "finalVoltage" and "finalAngle" are the calculated (through power flow or state estimation) voltage and angle of the bus. These parameters are optional key arguments that do not need to be given when the instance is created; in this case, they assume the 1 and 0 values ("flat start"). These values can be changed directly or through the runPowerFlow() method in the case class.
class bus:
	def __init__(self,number,name,PVtype,pLoad,qLoad,pGen,qGen,gsh,bsh,finalVoltage = None, finalAngle = None):
		self.number = int(number)
		self.name = str(name)
		self.PVtype = str(PVtype)
		self.pLoad = float(pLoad)
		self.qLoad = float(qLoad)
		self.pGen = float(pGen)
		self.qGen = float(qGen)
		self.bsh = float(bsh)
		self.gsh = float(gsh)

		self.finalVoltage = 1 if finalVoltage is None else finalVoltage
		self.finalAngle = 0 if finalAngle is None else finalAngle

	def __str__(self):
		tableformat = 'psql'
		print(' >> Bus	\'{0}\':'.format(self.name))
		tabRows = []
		tabHeader = ['Number', 'Name', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)', 'Active Generation\npGer (MW)', 'Reactive Generation\nqGen (MVAR)', 'Shunt capacitance\nbsh (p.u.)', 'Shunt conductance\ngsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (rad)']
		tabRows.append([self.number, self.name, self.PVtype, self.pLoad, self.qLoad, self.pGen, self.qGen, self.bsh, self.gsh, self.finalVoltage, self.finalAngle ])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

# (3) Branch object {{{1
# The branch object stores data for a particular branch of the system:
#--> "fromBus" and "toBus" are the numbers of respectively the first and second buses attached to the branch.
#--> "r"  and "x" are respectively the equivalent resistance and reactance of the branch.
#--> "bsh" is the shunt susceptance of the branch according to the pi model
#--> "a" is the turns ratio of the transformer attached to the branch when there is one.
class branch:
	def __init__(self,number,fromBus,toBus,r,x,bsh,a):
		self.number = int(number)
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
	def __init__(self,busNumber,ratedPower,ratedVoltage,H,D,ra,xL,xd,xPd,xPPd,tPdo,tPPdo,xq,xPq,xPPq,tPqo,tPPqo):
		self.busNumber = int(busNumber)
		self.ratedPower = float(ratedPower)
		self.ratedVoltage = float(ratedVoltage)
		self.H = float(H)
		self.D = float(D)
		self.ra = float(ra)
		self.xL = float(xL)
		self.xd = float(xd)
		self.xPd = float(xPd)
		self.xPPd = float(xPPd)
		self.tPdo = float(tPdo)
		self.tPPdo = float(tPPdo)
		self.xq = float(xq)
		self.xPq = float(xPq)
		self.xPPq = float(xPPq)
		self.tPqo = float(tPqo)
		self.tPPqo = float(tPPqo)

# (5) Fault object {{{1
class fault:
	def __init__(self,branch,location,openingTime):
		self.branch = branch
		self.location = location
		self.openingTime = openingTime
