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

import libraries.logmessages as logmgs

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
	def __init__(self, system_data, bus_data, branch_data, gen_data, machine_data = None, fault_data = None):
		self.name, self.Sb, self.Vb, self.droopType = system_data
		self.bus_data = bus_data
		self.branch_data = branch_data
		self.gen_data = gen_data
		self.machine_data = [] if machine_data is None else machine_data
		self.fault_data = [] if fault_data is None else fault_data

		# Obtaining the matrixes from the parameters:
		#	- These parameters are defined when the self.build functions are run
		#	- This means that the case class should be initialized AFTER the bus, branch and generator data are already loaded, that is, when busData, branchData and genData are already defined
		#	- This explains why in the load_case program the case class is only initialized after these data are defined and not before by using busData = [], branchData = [] and genData = []
		self.r, self.x, self.z, self.y, self.g, self.b, self.a, self.phi, self.gsh, self.bsh, self.K  = self.build_y()
		self.Y, self.G, self.B = self.build_Y();
		self.Ysh, self.Gsh, self.Bsh = self.build_Ysh();
		self.isP, self.isQ, self.isV, self.isT = self.build_is_matrices()

	def __str__(self):

		tableformat = 'psql' 

		print(logmgs.header_message('Printing case \'{0}\' data:'.format(self.name)))
		
		self.print_bus_data()
		self.print_branch_data()
		self.print_gen_data()
		self.print_fault_data()
		self.print_machine_data()
		return ''

	def get_nbus(self): return len(self.bus_data)
	def get_ngen(self): return len(self.gen_data)

	# runPowerFlow() runs the power flow of the case and updates the finalVoltage and finalAngle attributes of the bus instances in the busData list. If the power flow method is not successful it returns a false value.
	def run_power_flow(self):
		nbus = self.get_nbus()
		V, theta, r, itCount, success = pF.power_flow(self)
		if not success:
			print(' >>> Power flow update for case \'{0}\' instance returned not successful because the power flow method did not converge.'.format(self.name))
			V, theta = array(['NaN']*nbus), array(['NaN']*nbus)
			return False

		else:
			for k in range(nbus): self.busData[k].finalVoltage, self.busData[k].finalAngle = V[k], theta[k]

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
	def build_y(self):
		nbus = self.get_nbus()
		r = np.zeros((nbus,nbus))
		x = np.zeros((nbus,nbus))
		a = np.ones((nbus,nbus))
		gsh = np.zeros((nbus,nbus))
		bsh = np.zeros((nbus,nbus))
		phi = np.zeros((nbus,nbus))
		K = [ [] for i in range(nbus)]
		
		for branch in self.branch_data:
			load_bus_number = self.get_bus_number(branch.load_bus)
			tap_bus_number = self.get_bus_number(branch.tap_bus)

			K[load_bus_number] += [tap_bus_number]
			K[tap_bus_number] += [load_bus_number]

			r[load_bus_number, tap_bus_number] = branch.r
			r[tap_bus_number, load_bus_number] = branch.r

			x[load_bus_number, tap_bus_number] = branch.x
			x[tap_bus_number, load_bus_number] = branch.x

			gsh[load_bus_number, tap_bus_number] = branch.gsh	
			gsh[tap_bus_number, load_bus_number] = branch.gsh

			bsh[load_bus_number, tap_bus_number] = branch.bsh	
			bsh[tap_bus_number, load_bus_number] = branch.bsh

			phi[load_bus_number, tap_bus_number] = branch.phi

			a[load_bus_number, tap_bus_number] = branch.a

		for k in range(nbus): bsh[k,k] = self.bus_data[k].bsh
		
		z = array([[r[m,n] + 1j*x[m,n] for m in range(nbus)] for n in range(nbus)])
		y = array([[1/z[m,n] if z[m,n] != 0 else 0 for m in range(nbus)] for n in range(nbus)])
		g = np.real(copy(y))
		b = np.imag(copy(y))
		
		return r,x,z,y,g,b,a,phi,gsh,bsh,K

# Method case.buildY builds the admittance matrix Y of the system
	def build_Y(self):
		nbus = self.get_nbus()
		Y = np.zeros((nbus,nbus),dtype=complex)
		for k in range(nbus):
			for m in self.K[k]: Y[k,m] = -self.a[k,m]*np.exp(-1j*self.phi[k,m])*self.a[m,k]*np.exp(1j*self.phi[m,k])*(self.g[k,m] + 1j*self.b[k,m])
			Y[k,k] = self.bus_data[k].gsh + 1j*self.bus_data[k].bsh + sum([self.a[k,m]**2*(self.gsh[k,m] + 1j*self.bsh[k,m] + self.g[k,m] + 1j*self.b[k,m]) for m in self.K[k]])
		return Y, np.real(Y), np.imag(Y)

# Method case.buildY builds the shunt admittance matrix Y of the system
	def build_Ysh(self):
		nbus = self.get_nbus()
		Ysh = np.zeros((nbus,nbus),dtype=complex)
		for k in range(nbus):
			for m in self.K[k]: Ysh[k,m] = self.a[k,m]**2*(self.gsh[k,m] + 1j*self.bsh[k,m] + self.g[k,m] + 1j*self.b[k,m])
			Ysh[k,k] = sum(Ysh[k,m] for m in range(nbus))
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
		nbus = self.get_nbus()
		isP = np.ones((nbus,nbus))
		isQ = np.ones((nbus, nbus))
		isV = np.ones(nbus)
		isT = np.ones(nbus)

		for i in range(nbus):
			if self.bus_data[i].bus_type == 'VT':
				# In the case of VT buses, voltage magnitude and angles are fixed and known. The values used are the ones specified in the netfile as the bus final voltage and final angle values.
				isV[i] = 0
				isT[i] = 0
			if self.bus_data[i].bus_type == 'PV':
				isV[i] = 0
				isP[i,i] = 0
			if self.bus_data[i].bus_type in ['PQ', 'dPQ', 'PdQ']:
				isP[i,i] = 0
				isQ[i,i] = 0

		return isP, isQ, isV, isT
	
	def update_matrices(self):
			self.r, self.x, self.z, self.y, self.g, self.b, self.a, self.phi, self.gsh, self.bsh, self.K = self.build_y()
			self.Y, self.G, self.B = self.build_Y()
			self.Ysh, self.Gsh, self.Bsh = self.build_Ysh()
			self.isP, self.isQ, self.isV, self.isT = self.build_is_matrices()
		

# getBusNumber receives the name of a target bus and outputs its number in the busData list. If the name is not found in the list, the function returns an error and outputs the number -1
	def get_bus_number(self,target_bus_name):
		for bus in self.bus_data:
			if bus.name == target_bus_name: return bus.number
		else:
			raise NameError(' >> getBusNumber error: the provided target bus name was not found in the bus list.')
			return -1
			

# Function swapBuses is used to swap bus 'bus1' and bus 'bus2' in the busList and update the system matrixes and data lists to reflect that swapping.
	def swapBuses(self,bus1,bus2):
		i = self.get_bus_number(bus1)
		j = self.get_bus_number(bus2)
		temp_case = copy(self)
		temp_case.bus_data[i] = copy(self.bus_data[j])
		temp_case.bus_data[j] = copy(self.bus_data[i])
		
		# Resetting "from bus" and "to bus" numbers
		temp_case.bus_data[i].number = i
		temp_case.bus_data[j].number = j

		temp_case.update_matrixes()

		return temp_case

# The 'reduceMatrixes' method returns the matrixes of the reduced system. In order to do this, the system will need to be reorganized so that the buses attached to a generator are numbered first. Since this method alters the sequence of buses (that is, it reorganizes the number of the buses), it is recommended that this method be run through a copied instance, that is, run 'reducedCase = copy(case)' and then 'reducedCase.reduceGrid()'. Although this may seem like a problem, given that the user will not know which buses were re-numbered, the bus names defined won'be altered, meaning that the human-readable names do not change. This ultimately means that it is imperative to give each bus a human-readable name. 
# In order to do all this, the method creates a copy of the original case, called 'rCase' for 'reduced case', upon which it will operate and swap buses.
	def reduce_matrixes(self):
		# Building component matrixes
		ngen, nbus = self.get_ngen(), self.get_nbus()
		success = self.run_power_flow()
		#print('-'*50 + '\n REDUCING CASE {}\n'.format(self.name) + '-'*50 + '\n')

		# Checking if power flow is successful
		if not success: raise Exception(' Case \'{}\' matrix reduction not possible because the power flow calculations returned not successful'.format(self.name))

		# Creates a new case, rCase (short for "reduced case"). The idea is that all fixed power loads defined in the case are substituted for equivalent shunt admittances such that the initial conditions power flow is the same for both cases, but the reduced case can make use of the different 
		reduced_case = copy(self)
		reduced_case.name = 'Reduced ' + self.name

		# Moving generator buses to the top positions
		for i in range(nbus):
			for j in range(i+1, nbus):
				if reduced_case.is_gen(reduced.case.bus_data[j].name) and not reduced_case.is_gen(reduced_case.bus_data[i].name): reduced_case = reduced_case.swap_buses(reduced_case.bus_data[i].name, reduced_case.bus_data[j].name)

		# Calculating the equivalent conductance of the bus loads and adding those conductance to the bus shunt conductances
		YL = np.diag([ (reduced_case.bus_data[k].p_load  - 1j*reduced_case.bus_data[k].q_load)/reduced_case.bus_data[k].final_voltage**2 for k in range(reduced_case.get_nbus())])/reduced_case.Sb
		#print(YL)
		for k in range(nbus):
			reduced_case.bus_data[k].p_load, reduced_case.bus_data[k].q_load = 0, 0
			reduced_case.bus_data[k].gsh += np.real(YL[k, k])
			educed_casee.bus_data[k].bsh += np.imag(YL[k, k])

		reduced_case.update_matrixes()
		Y = copy(reduced_case.Y)

		Y1 = Y[0 : ngen, 0 : ngen]
		Y2 = Y[0 : ngen , ngen : nbus]
		Y3 = Y[ngen : nbus , 0 : ngen]
		Y4 = Y[ngen : nbus , ngen : nbus]

		Y_reduced = Y1 - Y2 @ inv(Y4) @ Y3

		# Reorganizing the first buses according to their generator order
		for i in range(ngen):
			for j in range(i, reduced_case.ngen):
				if reduced_case.get_bus_number(reduced_case.gen_data[i].bus_name) > reduced_case.get_bus_number(reduced_case.gen_data[j].bus_name) : reduced_case.gen_data[i], reduced_case.gen_data[j] = reduced_case.gen_data[j], reduced_case.gen_data[i]
		return Y_reduced, r_case

	def print_bus_data(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' bus list'.format(self.name)))
		tab_rows = []
		tab_header = ['Number', 'Name', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)', 'Shunt conductance\ngsh (p.u.)', 'Shunt susceptancee\nbsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (deg)']
		tab_header = [logmgs.bold_string(x) for x in tab_header]
		for bus in self.bus_data: tab_rows.append([bus.number, bus.name, bus.bus_type, bus.p_load, bus.q_load, bus.gsh, bus.bsh, bus.final_voltage, bus.final_angle])

		print(tabulate(tab_rows,headers=tab_header, numalign='right', tablefmt=tableformat))

	def print_branch_data(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' branch list'.format(self.name)))
		tab_rows = []
		tab_header = ['Number','Load Bus\n(\"From\" bus)', 'Tap Bus\n(\"To\" bus)', 'Resistance\n r (p.u.)', 'Reactance\n x (p.u.)', 'Shunt conductance\n gsh (p.u.)', 'Shunt susceptance\n bsh (p.u.)', 'Transformer\nturns ratio (a)', 'Transformer\nphase shift (phi)', 'Active\nPower Transfer (MW)', 'Reactive\nPower Transfer (MVAR)', 'Active\nPower Loss (MW)', 'Reactive\nPower Loss (MVAR)', 'Reactive\nShunt Generated\nPower (MVAR)']
		for branch in self.branch_data: tab_rows.append([branch.number, branch.load_bus, branch.tap_bus, branch.r, branch.x, branch.gsh, branch.bsh, branch.a,branch.phi*180/np.pi, branch.active_transfer*self.Sb, branch.reactive_transfer*self.Sb, branch.active_loss*self.Sb, branch.reactive_loss*self.Sb, branch.shunt_reactive_power*self.Sb])

		print(tabulate(tab_rows,headers=tab_header, numalign='right', tablefmt=tableformat))

	def print_gen_data(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql' 

		print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' generator list (sorted by bus number)'.format(self.name)))
		tab_rows = []
		tab_header = ['Bus', 'Rated Active \n Power (MW)', 'Rated Reactive \n Power(MVAr)', 'Rated Voltage \n (kV)', 'Active Power\nGenerated (MW)', 'Reactive Power\nGenerated (MVAr)', 'Minimum Active\nGeneration (MW)', 'Maximum Active\n Generation (MVAr)', 'Minimum Reactive\nGeneration (MVAr)', 'Maximum Reactive\n Generation (MVAr)', 'Reference Voltage\n(kV)', 'Active Droop\nRamp (MW/kV)', 'Reactive Droop\nRamp (MVAr/kV)']
		for gen in self.gen_data:
			tab_rows.append([gen.bus_name, gen.p_rated, gen.q_rated, gen.v_rated, gen.p_gen, gen.q_gen, gen.p_min, gen.p_max, gen.q_min, gen.q_max, gen.v_ref, gen.p_ramp, gen.q_ramp])

		print(tabulate(tab_rows,headers=tab_header, numalign='right', tablefmt=tableformat))

	def print_fault_data(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql'
 
		if self.fault_data == [] : print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' has no fault data.'.format(self.name)))
		else:
			print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\'  possible faults list'.format(self.name)))
			tab_rows = []
			tab_header = ['Branch Number', 'Location', 'Opening Time']
			for fault in self.fault_data: tab_rows.append([fault.branch, fault.location, fault.opening_time])
			print(tabulate(tab_rows,headers=tab_header, numalign='right', tablefmt=tableformat))

	def print_machine_data(self,**kwargs):
		if 'tablefmt' in kwargs: tableformat = kwargs['tablefmt']
		else: tableformat = 'psql'
 
		if self.machine_data == [] : print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' has no machine data.'.format(self.name)))
		else:
			print('\n' + logmgs.green_indicator_string(' >>', ' Case \'{0}\' machine data'.format(self.name)))
			tab_rows = []
			tab_header = ['Bus']
			for machine in self.machine_data: tab_rows.append([machine.bus])
			print(tabulate(tab_rows,headers=tab_header, numalign='right', tablefmt=tableformat))

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
	def __init__(self, number, name, bustype, final_voltage, final_angle, pLoad, qLoad, gsh, bsh, vmax, vmin):
		self.number = int(number)
		self.name = str(name)
		self.bus_type = str(bustype)
		self.final_voltage = 1 if final_voltage is None else final_voltage
		self.final_angle = 0 if final_angle is None else final_angle
		self.p_load = float(pLoad)
		self.q_load = float(qLoad)
		self.bsh = float(bsh)
		self.gsh = float(gsh)
		self.v_max = float(vmax)
		self.v_min = float(vmin)
		self.v_flag = False

	def __str__(self):
		tableformat = 'psql'
		print(' >> Bus	\'{0}\':'.format(self.name))
		tabRows = []
		tabHeader = ['Number', 'Name', 'Type', 'Active Load\npLoad (MW)', 'Reactive Load\nqLoad (MVAR)','Shunt capacitance\nbsh (p.u.)', 'Shunt conductance\ngsh (p.u.)', 'Final voltage\nV (p.u.)', 'Final angle\ntheta (deg)', 'Maximum\nVoltage\n(p.u.)', 'Minimum\nVoltage\n(p.u.)']
		tabRows.append([self.number, self.name, self.bus_type, self.pLoad, self.qLoad, self.bsh, self.gsh, self.finalVoltage, self.finalAngle*180/np.pi, self.vmax, self.vmin ])
		return tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat)

# Function get_bus takes a bus name and returns the bus object that has that name.
def get_bus(bus_name, bus_list):
	for bus in bus_list:
		if bus_name == bus.name: return bus
	else: return None
		
# (3) Branch object {{{1
# The branch object stores data for a particular branch of the system:
#--> "fromBus" and "toBus" are the numbers of respectively the first and second buses attached to the branch.
#--> "r"  and "x" are respectively the equivalent resistance and reactance of the branch.
#--> "bsh" is the shunt susceptance of the branch according to the pi model
#--> "a" is the turns ratio of the transformer attached to the branch when there is one.
#--> "phi" is the a
class branch:
	def __init__(self, number, load_bus, tap_bus, r, x, gsh, bsh, a, phi, minphi, maxphi, status, activeTransfer = None, reactiveTransfer = None, activeLoss = None, reactiveLoss = None, shuntReactivePower = None):
		self.number = int(number)
		self.load_bus = str(load_bus)
		self.tap_bus = str(tap_bus)
		self.r = float(r)
		self.x = float(x)
		self.gsh = float(gsh)
		self.bsh = float(bsh)
		self.a = float(a)
		self.phi = float(phi)
		self.minphi = float(minphi)
		self.maxphi = float(maxphi)
		self.status = bool(status)
		self.reactive_transfer = float('NaN') if reactiveTransfer is None else float(reactiveTransfer)
		self.active_transfer = float('NaN') if activeTransfer is None else float(activeTransfer)
		self.active_loss = float('NaN') if activeLoss is None else float(activeLoss)
		self.reactive_loss = float('NaN') if reactiveLoss is None else float(reactiveLoss)
		self.shunt_reactive_power = float('NaN') if shuntReactivePower is None else float(shuntReactivePower)

	def __str__(self):
		tableformat = 'psql'
		print(' >> Branch from bus \'{0}\' to bus \'{1}\':'.format(self.fromBus, self.toBus))
		tabRows = []
		tabHeader = ['Number', 'Tap bus', 'Load bus', 'Resistance\n(p.u.)', 'Reactance\n(p.u.)','Shunt conductance\ngsh (p.u.)', 'Shunt susceptance\nbsh (p.u.)', 'Transformer\nRatio (-)', 'Transformer phase\nangle (deg)', 'Minimum phase\nangle (p.u.)', 'Maximum phase\nangle (p.u.)', 'Status', 'Active Power\nTransfer (p.u.)', 'Reactive Power\nTransfer (p.u.)', 'Active Power\nLoss (p.u.)', 'Reactive Power\nLoss (p.u.)', 'Shunt reactive\nload (MVAr)']
		tabRows.append([self.number, self.toBus, self.fromBus, self.r, self.x, self.gsh, self.bsh, self.a, self.phi*180/np.pi, self.minphi*180/np.pi, self.maxphi*180/pi, self.status, self.activeTransfer, self.reactiveTransfer, self.activeLoss, self.reactiveLoss, self.shuntReactivePower])
#
		return tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat)


# (4) Generator object {{{1
class generator:
	def __init__(self, bus_name, p_gen, q_gen, q_max, q_min, p_max, p_min, v_ref, p_rated, q_rated, v_rated, p_ramp, q_ramp):#,H,D,ra,xL,xd,xPd,xPPd,tPdo,tPPdo,xq,xPq,xPPq,tPqo,tPPqo, Ke, Te, vRef, KPss, Tw, T1, T2, tG, tT, modelDepth):
		self.bus_name = str(bus_name)
		self.p_rated = str(p_rated)
		self.q_rated = str(q_rated)
		self.v_rated = str(v_rated)
		self.p_gen = str(p_gen)
		self.q_gen = str(q_gen)
		self.p_min = str(p_min)
		self.p_max = str(p_max)
		self.q_min = str(q_min)
		self.q_max = str(q_max)
		self.v_ref = str(v_ref)
		self.p_ramp = str(p_ramp)
		self.q_ramp = str(q_ramp)

	def __str__(self):
		tableformat = 'psql' 

		print('\n >> Generator at bus \'{}\''.format(self.bus_name))
		tabRows = []
		tabHeader = ['Generated P (MW)', 'Generated Q (MVAr)', 'Max Q (MVAr)', 'Min Q (MVAr)', 'Max P (MW)', 'Min P (MW)', 'Reference V (p.u.)', 'Rated P (MW)', 'Rated Q (MVAr)', 'Rated V (kV)', 'P ramp (MW/kV)', 'Q ramp (MVAr/kV)']
		gen = self
		tabRows.append([self.p_gen, self.q_gen, self.q_max, self.q_min, self.p_max, self.p_min, self.v_ref, self.p_rated, self.q_rated, self.v_rated, self.p_ramp, self.q_ramp])

		print(tabulate(tabRows,headers=tabHeader, numalign='right', tablefmt=tableformat))
		return ''

# is_gen takes a string bus_name and a list of generator instances; it sweeps the list to check if, in the provided generator list, there is a generator attached to the bus. If so, it returns the generator instance; if not, it returns false.
def is_gen(bus_name, gen_list):
	for gen in gen_list:
		if bus_name == gen.bus_name: return gen
	return False



# (5) Machine object {{{1
class machine:
	def __init__(self):
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
