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

# (1) Case object {{{1
# The case object stores data about a particular "case". "Cases" are a strain of data which comprise simulation parameters of a given system. These parameters include topological data like bus, branch and generator data and simulatory parameters, such as simulation time, plot data, tolerances et cetera. This object serves the purpose of giving the user an easy way to modify and customize the set of parameters -- be them topological or simulatory -- without having to generate a whole new file. In general, a case file is built from a *.net supplied by the user; see the loadCase function in file loadCase.py
# --> "name" is a human-readable name that can be remembered by the user;
# --> "busData" is a list of "bus" objects (see (2)), which store information about buses like injected power, load power, name, number;
# --> "branchData" is a list of "branch" objects (see (3)) which store information about the branches of the system, such as resistance and reactance, transformer turns ratio and so on;
# --> "genData" is a list of "generator" objects (see (4)) which sotre information about system generators. These generators are modelled as synchronous machines; the "generator" object sotres basically constructive parameters of the machines.
# --> "faultData" is a list of "fault" objects (see (5)) which store information about possible faults in the system. The listed faults will be used for dynamical simulation.
class case:
	def __init__(self,name, busData,branchData,genData,faultData,Sb,Vb):
		self.name = name
		self.busData = busData
		self.branchData = branchData
		self.genData = genData
		self.faultData = faultData	
		self.Sb = Sb			# Base power value
		self.Vb = Vb			# Base voltage value
		self.y = 0
		self.Y = 0

	def buildMatrixes(self):
		nBus = len(self.branchData)
		z = array([[r[m,n] + 1j*x[m,n] for m in range(nBus)] for n in range(nBus)])
		y = array([[1/z[m,n] if z[m,n] != 0 else 0 for m in range(nBus)] for n in range(nBus)])

		g = real(copy(y))
		b = imag(copy(y))

		bsh = array([[ bsh[m,n]/2 if m!=n else bsh[m,n] for n in range(nBus)] for m in range(nBus)])

		Y  = -a*transpose(a)*y

		for k in range(nBus): Y[k,k] = 1j*bsh[k,k] + sum([a[k,m]**2*y[k,m] + 1j*bsh[k,m] for m in K[k]])

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
	def __init__(self,fromBus,toBus,neighbors,r,x,bsh,a):
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

# (5) Fault object {{{1
class fault:
	def __init__(self,branch,location,openingTime):
		self.branch = branch
		self.location = location
		self.openingTime = openingTime
