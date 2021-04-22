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
# --> Y is the multi-dimensional thevenin equivalent conductance of the grid.
# --> Y is the multi-dimensional thevenin equivalent of the shunt conductances of the grid.

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

		self.r, self.x, self.z, self.y, self.g, self.b, self.a, self.phi, self.gsh, self.bsh, self.K  = self.buildy()		
		#self.z = array([[self.r[m,n] + 1j*self.x[m,n] for m in range(self.nBus)] for n in range(self.nBus)])
		#self.y = array([[1/self.z[m,n] if self.z[m,n] != 0 else 0 for m in range(self.nBus)] for n in range(self.nBus)])
		#self.g = real(copy(self.y))
		#self.b = imag(copy(self.y))
		#self.bsh = array([[ self.bsh[m,n]/2 if m!=n else self.bsh[m,n] for n in range(self.nBus)] for m in range(self.nBus)])

		self.Y, self.G, self.B = self.buildY();
		self.Ysh, self.Gsh, self.Bsh = self.buildYsh();
		self.isP, self.isQ, self.isV, self.isT = self.build_is_matrices()

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
		V, theta, r, itCount, success = pF.powerFlow(self)
		if not success:
			print(' >>> Power flow update for case class instance returned not successful because the power flow method did not converge.')
			V, theta = array(['None']*self.nBus), array(['None']*self.nBus)
			return False

		else:
			for k in range(self.nBus): self.busData[k].finalVoltage, self.busData[k].finalAngle = V[k], theta[k]

			for k in range(self.nBus):
				self.busData[k].pGen  = V[k]**2*self.G[k,k] + V[k]*sum([ V[m]*(self.G[k,m]*cos(theta[k] - theta[m]) + self.B[k,m]*sin(theta[k] - theta[m])) for m in self.K[k]]) + self.busData[k].pLoad/self.Sb
				self.busData[k].pGen = self.Sb*np.round(self.busData[k].pGen*1e10)/1e10
				self.busData[k].qGen = -V[k]**2*self.B[k,k] + V[k]*sum([ V[m]*(self.G[k,m]*sin(theta[k] - theta[m]) - self.B[k,m]*cos(theta[k] - theta[m])) for m in self.K[k]]) + self.busData[k].qLoad/self.Sb
				self.busData[k].qGen = self.Sb*np.round(self.busData[k].qGen*1e10)/1e10
				self.busData[k].finalAngle = 180/pi*theta[k]
			
			for k in range(len(self.branchData)):
				j = self.getBusNumber(self.branchData[k].fromBus)
				n = self.getBusNumber(self.branchData[k].toBus)
				self.branchData[k].activeTransfer = (self.a[j,n]*V[j])**2*self.g[j,n] - self.a[j,n]*self.a[n,j]*V[j]*V[n]*(self.g[j,n]*cos(theta[j] - theta[n] + self.phi[j,n] - self.phi[n,j]) + self.b[j,n]*sin(theta[j] - theta[n] + self.phi[j,n] - self.phi[n,j]))
				self.branchData[k].reactiveTransfer = -(self.a[j,n]*V[j])**2*(self.b[j,n] + self.bsh[j,n]) + self.a[j,n]*self.a[n,j]*V[j]*V[n]*( -self.g[j,n]*sin(theta[j] - theta[n] + self.phi[j,n] - self.phi[n,j]) + self.b[j,n]*cos(theta[j] - theta[n] + self.phi[j,n] - self.phi[n,j]) )
				self.branchData[k].activeLoss = self.g[j,n]*np.abs(V[j]*np.exp(1j*theta[j]) - V[n]*np.exp(1j*theta[n]))**2
				self.branchData[k].reactiveLoss = -self.b[j,n]*np.abs(V[j]*np.exp(1j*theta[j]) - V[n]*np.exp(1j*theta[n]))**2
				self.branchData[k].shuntReactivePower = -self.bsh[j,n]*(V[j]**2 + V[n]**2)

			return True


# Function case.update_matrixes is meant to update matrixes r,x,a,bsh and K everytime these are called, or everytime they are needed. This is done to prevent inconsistencies if the user changes a variable directly, say for example:
# case.branchData[1].r = 5
# In this case the user has directly updated the value of a branch resistance; all matrixes must be re-calculated. Function update_matrixes is called whenever case.r, case.x, case.a, case.bsh, case.K are called, so they are recalculated for the present values of bus and branch data.
	def buildy(self):
		r = np.zeros((self.nBus,self.nBus))
		x = np.zeros((self.nBus,self.nBus))
		a = np.ones((self.nBus,self.nBus))
		gsh = np.zeros((self.nBus,self.nBus))
		bsh = np.zeros((self.nBus,self.nBus))
		phi = np.zeros((self.nBus,self.nBus))
		K = [ [] for i in range(self.nBus)]
		
		for branch in self.branchData:
			fromBusN = self.getBusNumber(branch.fromBus)
			toBusN = self.getBusNumber(branch.toBus)

			K[fromBusN] += [toBusN]
			K[toBusN] += [fromBusN]

			r[fromBusN, toBusN] = branch.r
			r[toBusN, fromBusN] = branch.r

			x[fromBusN, toBusN] = branch.x
			x[toBusN, fromBusN] = branch.x

			gsh[fromBusN, toBusN] = branch.gsh	
			gsh[toBusN, fromBusN] = branch.gsh

			bsh[fromBusN, toBusN] = branch.bsh	
			bsh[toBusN, fromBusN] = branch.bsh

			phi[fromBusN, toBusN] = branch.phi

			a[fromBusN, toBusN] = branch.a

		for k in range(self.nBus): bsh[k,k] = self.busData[k].bsh
		
		z = array([[r[m,n] + 1j*x[m,n] for m in range(self.nBus)] for n in range(self.nBus)])
		y = array([[1/z[m,n] if z[m,n] != 0 else 0 for m in range(self.nBus)] for n in range(self.nBus)])
		g = np.real(copy(y))
		b = np.imag(copy(y))
		
		return r,x,z,y,g,b,a,phi,gsh,bsh,K

# Method case.buildY builds the admittance matrix Y of the system
	def buildY(self):
		Y = np.zeros((self.nBus,self.nBus),dtype=complex)
		for k in range(self.nBus):
			for m in self.K[k]: Y[k,m] = -self.a[k,m]*np.exp(-1j*self.phi[k,m])*self.a[m,k]*np.exp(1j*self.phi[m,k])*(self.g[k,m] + 1j*self.b[k,m])
			Y[k,k] = self.busData[k].gsh + 1j*self.busData[k].bsh + sum([self.a[k,m]**2*(self.gsh[k,m] + 1j*self.bsh[k,m] + self.g[k,m] + 1j*self.b[k,m]) for m in self.K[k]])
		return Y, np.real(Y), np.imag(Y)

# Method case.buildY builds the shunt admittance matrix Y of the system
	def buildYsh(self):
		Ysh = np.zeros((self.nBus,self.nBus),dtype=complex)
		for k in range(self.nBus):
			for m in self.K[k]: Ysh[k,m] = self.a[k,m]**2*(self.gsh[k,m] + 1j*self.bsh[k,m] + self.g[k,m] + 1j*self.b[k,m])
			Ysh[k,k] = sum(Ysh[k,m] for m in range(self.nBus))
		return Ysh, np.real(Ysh), np.imag(Ysh)

# Method build_is builds the isP, isQ, isV and isT arrays
	def build_is_matrices(self):
		# isP, isQ and isV are the matrixes/array that flag P, Q and V measures.
		# isP[n,n] == 0 means that the n-th bus active power is not a variable for the power flow method, that is, has a fixed power injection. This generally happens for PQ and PV bus types. The same goes for isQ[n,n] == 0.
		# isP[n,n] == 1 means that the n-th bus active power is a variable function of V and theta, and should be estimated. This is generally the case for VT bus types. The same goes for isQ[n,n].
		# isV[n] == 0 means that the voltage of the n-th bus is not a variable and should not be calculated, that is, has a fixed value. This is generally the case for PQ bus types.
		# isV[n] == 1 means that the voltage of the n-th bus is a varuable and should be estimated. This is generally the case for PQ buses.
		# isT[n] == 0 means that the voltage angle of the n-th bus is not a variable and should not be estimated, that is, has a fixed value. This generally happens for VT buses.
		# isT[n] == 1 means that the angle of the n-th bus is not a variable and should be estimated. This is generally true for PQ and PV buses.
		isP = np.eye(self.nBus)
		isQ = np.eye(self.nBus)
		isV = np.ones(self.nBus)
		isT = np.ones(self.nBus)
		
		for i in range(self.nBus):
			if self.busData[i].PVtype == 'VT':
				# In the case of VT buses, voltage magnitude and angles are fixed and known. The values used are the ones specified in the netfile as the bus final voltage and final angle values.
				isV[i] = 0
				isT[i] = 0
				isP[i,i] = 0
				isQ[i,i] = 0
			if self.busData[i].PVtype == 'PV':
				isV[i] = 0
				isP[i,i] = 1
				isQ[i,i] = 0
			if self.busData[i].PVtype == 'PQ':
				isP[i,i] = 1
				isQ[i,i] = 1

		return isP, isQ, isV, isT
	
	def update_matrixes(self):
			self.r, self.x, self.z, self.y, self.g, self.b, self.a, self.phi, self.gsh, self.bsh, self.K = self.buildy()
			self.Y, self.G, self.B = self.buildY()
			self.Ysh, self.Gsh, self.Bsh = self.buildYsh()
			self.isP, self.isQ, self.isV, self.isT = self.build_is_matrices()
		

# getBusNumber receives the name of a target bus and outputs its number in the busData list. If the name is not found in the list, the function returns an error and outputs the number -1
	def getBusNumber(self,targetBusName):
		for bus in self.busData:
			if bus.name == targetBusName: return bus.number
		else:
			raise NameError(' >> getBusNumber error: the provided target bus name was not found in the bus list.')
			return -1
			

# Function swapBuses is used to swap bus 'bus1' and bus 'bus2' in the busList and update the system matrixes and data lists to reflect that swapping.
	def swapBuses(self,bus1,bus2):
		i = self.getBusNumber(bus1)
		j = self.getBusNumber(bus2)
		tempCase = copy(self)
		tempCase.busData[i] = copy(self.busData[j])
		tempCase.busData[j] = copy(self.busData[i])
		
		# Resetting "from bus" and "to bus" numbers
		tempCase.busData[i].number = i
		tempCase.busData[j].number = j

		tempCase.update_matrixes()

		return tempCase

# isGen(busN) returns True if the bus with number busN has a generator attached to it.
	def isGen(self, busName):
		for gen in self.genData:
			if gen.busName == busName: return True

		return False

# The 'reduceMatrixes' method returns the matrixes of the reduced system. In order to do this, the system will need to be reorganized so that the buses attached to a generator are numbered first. Since this method alters the sequence of buses (that is, it reorganizes the number of the buses), it is recommended that this method be run through a copied instance, that is, run 'reducedCase = copy(case)' and then 'reducedCase.reduceGrid()'. Although this may seem like a problem, given that the user will not know which buses were re-numbered, the bus names defined won'be altered, meaning that the human-readable names do not change. This ultimately means that it is imperative to give each bus a human-readable name. 
# In order to do all this, the method creates a copy of the original case, called 'rCase' for 'reduced case', upon which it will operate and swap buses.
	def reduceMatrixes(self):
		# Building component matrixes
		nGen, nBus = self.nGen, self.nBus
		success = self.runPowerFlow()
		#print('-'*50 + '\n REDUCING CASE {}\n'.format(self.name) + '-'*50 + '\n')

		# Checking if power flow is successful
		if not success: raise Exception(' Case \'{}\' matrix reduction not possible because the power flow calculations returned not successful'.format(self.name))

		# Creates a new case, rCase (short for "reduced case"). The idea is that all fixed power loads defined in the case are substituted for equivalent shunt admittances such that the initial conditions power flow is the same for both cases, but the reduced case can make use of the different 
		rCase = copy(self)
		rCase.name = 'Reduced ' + self.name

		# Moving generator buses to the top positions
		for i in range(nBus):
			for j in range(i+1, nBus):
				if rCase.isGen(rCase.busData[j].name) and not rCase.isGen(rCase.busData[i].name):
					rCase = rCase.swapBuses(rCase.busData[i].name, rCase.busData[j].name)

		# Calculating the equivalent conductance of the bus loads and adding those conductance to the bus shunt conductances
		YL = np.diag([ (rCase.busData[k].pLoad  - 1j*rCase.busData[k].qLoad)/rCase.busData[k].finalVoltage**2 for k in range(rCase.nBus)])/rCase.Sb
		#print(YL)
		for k in range(rCase.nBus):
			rCase.busData[k].pLoad, rCase.busData[k].qLoad = 0, 0
			rCase.busData[k].gsh += np.real(YL[k, k])
			rCase.busData[k].bsh += np.imag(YL[k, k])

		rCase.update_matrixes()
		Y = copy(rCase.Y)

		Y1 = Y[0 : nGen, 0 : nGen]
		Y2 = Y[0 : nGen , nGen : nBus]
		Y3 = Y[nGen : nBus , 0 : nGen]
		Y4 = Y[nGen : nBus , nGen : nBus]

		Yred = Y1 - Y2 @ inv(Y4) @ Y3

		# Reorganizing the first buses according to their generator order
		for i in range(rCase.nGen):
			for j in range(i, rCase.nGen):
				if rCase.getBusNumber(rCase.genData[i].busName) > rCase.getBusNumber(rCase.genData[j].busName) : rCase.genData[i], rCase.genData[j] = rCase.genData[j], rCase.genData[i]
		return Yred, rCase

	def printBusData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' bus list'.format(self.name))
		tabRows = []
		tabHeader = ['Number', 'Name', 'Load Area', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)', 'Active Generation\npGer (MW)', 'Reactive Generation\nqGen (MVAR)', 'Shunt conductance\ngsh (p.u.)', 'Shunt susceptancee\nbsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (deg)']
		for bus in self.busData: tabRows.append([bus.number, bus.name, bus.loadArea, bus.PVtype, bus.pLoad, bus.qLoad, bus.pGen, bus.qGen, bus.gsh, bus.bsh, bus.finalVoltage, bus.finalAngle])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

	def printBranchData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' branch list'.format(self.name))
		tabRows = []
		tabHeader = ['Number','From Bus\n(Tap bus)', 'To Bus\n(Z bus)', 'Resistance\n r (p.u.)', 'Reactance\n x (p.u.)', 'Shunt conductance\n gsh (p.u.)', 'Shunt susceptance\n bsh (p.u.)', 'Transformer\nturns ratio (a)', 'Transformer\nphase shift (phi)', 'Active\nPower Transfer (MW)', 'Reactive\nPower Transfer (MVAR)', 'Active\nPower Loss (MW)', 'Reactive\nPower Loss (MVAR)', 'Reactive\nShunt Generated\nPower (MVAR)']
		for branch in self.branchData: tabRows.append([branch.number, branch.fromBus, branch.toBus, branch.r, branch.x, branch.gsh, branch.bsh, branch.a,branch.phi*180/np.pi, branch.activeTransfer*self.Sb, branch.reactiveTransfer*self.Sb, branch.activeLoss*self.Sb, branch.reactiveLoss*self.Sb, branch.shuntReactivePower*self.Sb])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))

	def printGenData(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n >> Case \'{0}\' generator list (sorted by bus number)'.format(self.name))
		tabRows = []
		tabHeader = ['Bus', 'Rated Power (MW)', 'H (p.u.)', 'D (p.u.)', 'ra (p.u.)', 'xL (p.u.)', 'xd (p.u.)', 'xPd (p.u.)', 'xPPd (p.u.)', 'tPdo (s)', 'tPPdo (s)', 'xq (p.u.)', 'xPq (p.u.)', 'xPPq (p.u.)', 'tPqo (s)', 'tPPqo (s)', 'Ke (-)', 'Te (s)', 'vRef (p.u.)', ' KPss (-)', 'T1 (s)', 'T1 (s)', 'tG (s)', 'tT (s)', 'Depth']
		for gen in self.genData:
			tabRows.append([gen.busName, gen.ratedPower,gen.H, gen.D, gen.ra, gen.xL, gen.xd, gen.xPd, gen.xPPd, gen.tPdo, gen.tPPdo, gen.xq, gen.xPq, gen.xPPq, gen.tPqo, gen.tPPqo, gen.Ke, gen.Te, gen.vRef, gen.KPss, gen.T1, gen.T2, gen.tG, gen.tT, gen.modelDepth])

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
	def __init__(self,number,name,loadArea,PVtype,pLoad,qLoad,pGen,qGen,gsh,bsh,finalVoltage = None, finalAngle = None):
		self.number = int(number)
		self.name = str(name)
		self.PVtype = str(PVtype)
		self.pLoad = float(pLoad)
		self.qLoad = float(qLoad)
		self.pGen = float(pGen)
		self.qGen = float(qGen)
		self.bsh = float(bsh)
		self.gsh = float(gsh)

		self.loadArea = int(loadArea)

		self.finalVoltage = 1 if finalVoltage is None else finalVoltage
		self.finalAngle = 0 if finalAngle is None else finalAngle

	def __str__(self):
		tableformat = 'psql'
		print(' >> Bus	\'{0}\':'.format(self.name))
		tabRows = []
		tabHeader = ['Number', 'Name', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)', 'Active Generation\npGer (MW)', 'Reactive Generation\nqGen (MVAR)', 'Shunt capacitance\nbsh (p.u.)', 'Shunt conductance\ngsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (rad)']
		tabRows.append([self.number, self.name, self.PVtype, self.pLoad, self.qLoad, self.pGen, self.qGen, self.bsh, self.gsh, self.finalVoltage, self.finalAngle ])
#
		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))
#
# (3) Branch object {{{1
# The branch object stores data for a particular branch of the system:
#--> "fromBus" and "toBus" are the numbers of respectively the first and second buses attached to the branch.
#--> "r"  and "x" are respectively the equivalent resistance and reactance of the branch.
#--> "bsh" is the shunt susceptance of the branch according to the pi model
#--> "a" is the turns ratio of the transformer attached to the branch when there is one.
#--> "phi" is the a
class branch:
	def __init__(self,number,fromBus,toBus,r,x,gsh,bsh,a,phi,activeTransfer = 0,reactiveTransfer = 0, activeLoss = 0, reactiveLoss = 0, shuntReactivePower = 0):
		self.number = int(number)
		self.fromBus = str(fromBus)
		self.toBus = str(toBus)
		self.r = float(r)
		self.x = float(x)
		self.gsh = float(gsh)
		self.bsh = float(bsh)
		self.a = float(a)
		self.phi = float(phi)
		self.reactiveTransfer = float(reactiveTransfer)
		self.activeTransfer = float(activeTransfer)
		self.activeLoss = float(activeLoss)
		self.reactiveLoss = float(reactiveLoss)
		self.shuntReactivePower = float(shuntReactivePower)

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
	def __init__(self,busName,ratedPower,ratedVoltage,H,D,ra,xL,xd,xPd,xPPd,tPdo,tPPdo,xq,xPq,xPPq,tPqo,tPPqo, Ke, Te, vRef, KPss, Tw, T1, T2, tG, tT, modelDepth):
		self.busName = str(busName)
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

		self.pm0 = 0
		self.efd0 = 0
		self.el0 = 0
		self.delta0 = 0
		self.omega0 = 0
		self.vAVR0 = 0
		self.vPSS0 = 0
		self.vWash0 = 0
		self.pPm0 = 0	# Initial value of the derivative of the mechanical power
		self.Ke = float(Ke)
		self.Te = float(Te)
		self.KPss = float(KPss)
		self.Tw = float(Tw)
		self.T1 = float(T1)
		self.T2 = float(T2)
		self.tG = float(tG)
		self.tT = float(tT)
		self.vRef = float(vRef)
		self.modelDepth = int(modelDepth)

		self.kQ = 1
		self.kP = 1
		self.P0	= 1
		self.Q0 = 0
		self.w0 = 0
		self.V0 = 0
		self.kReg = 0.01

	def __str__(self):
		tableformat = 'psql' 

		print('\n >> Generator at bus number {}'.format(self.busNumber))
		tabRows = []
		tabHeader = ['Bus', 'Rated Power (MW)', 'H (p.u.)', 'D (p.u.)', 'ra (p.u.)', 'xL (p.u.)', 'xd (p.u.)', 'xPd (p.u.)', 'xPPd (p.u.)', 'tPdo (s)', 'tPPdo (s)', 'xq (p.u.)', 'xPq (p.u.)', 'xPPq (p.u.)', 'tPqo (s)', 'tPPqo (s)', 'Ke (-)', 'Te (s)', 'vRef (p.u.)', ' KPss (-)', 'T1 (s)', 'T1 (s)', 'tG (s)', 'tT (s)', 'Depth']
		gen = self
		tabRows.append([ self.busNumber, gen.ratedPower,gen.H, gen.D, gen.ra, gen.xL, gen.xd, gen.xPd, gen.xPPd, gen.tPdo, gen.tPPdo, gen.xq, gen.xPq, gen.xPPq, gen.tPqo, gen.tPPqo, gen.Ke, gen.Te, gen.vRef, gen.KPss, gen.T1, gen.T2, gen.tG, gen.tT, gen.modelDepth])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))
		return ''


# (5) Fault object {{{1
class fault:
	def __init__(self,branch,location,openingTime):
		self.branch = branch
		self.location = location
		self.openingTime = openingTime
