# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: loadCase function
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.2.1
# DESCRIPTION: # loadcase is a function that takes the pointer "fileName" to a file and returns a case structure (see (2.6)).
# The structure is composed of bus, branch, generator and fault data, comprised of the lists of
# buses, branches, generators and faults. Each of these have a particular set of attributes that
# are explained in their respective object definitions.
# -------------------------------------------------

import libraries.logmessages as logmgs
import numpy as np
import libraries.classes as cL

from copy import deepcopy
copy = deepcopy

# PRELIMINARY FUNCTIONS ----------------------------------------------------------------------- {{{1
def isfloat(x):
	try: float(x)
	except ValueError: return False
	return True

def isinteger(x):
	try: int(x)
	except ValueError: return False
	return True

# VariableCheckError is used when a given variable, as declared in the netfile, is not in the type or range it should be. Additionally, it is also used to confirm incosistencies between bus and generator types. {{{2
class VariableCheckError(Exception):
	def __init__(self, *args):
		self.message = args[0]
	def __str__(self):
		return logmgs.error_message(self.message)

# Definitions of NaN, posinf and neginf, true, false and bool: these return arrays of possible definitions of NaN (Not Available Now or Not A Number), positive infinity and negative infinity, boolean true or false and boolean. {{{2
# These functions work as a roster of possible definitions for these terms that the user can input in the netfile, and are primarily used in the variable_check function to check if a variable is defined as one of these words.
# All definitions automatically become protected keywords that have to be expressly allowed in order for the check to pass. Hence, adding keywords to these rosters should be done with some care.
def nan_definitions(): return ['NaN', 'nan']
def posinf_definitions(): return ['inf', 'posinf', '+inf', 'infinite', 'infinity', 'INF']
def neginf_definitions(): return ['-inf', 'neginf', '-infinite', '-infinity', '-INF']
def true_definitions(): return ['1', 'true', 'TRUE', 'on', 'ON'] 
def false_definitions(): return ['0', 'false', 'FALSE', 'off', 'OFF']
def bool_definitions(): return true_definitions() + false_definitions()

# isbool checks if a string is a boolean according to bool_definitions() {{{2
def isbool(x): return x in bool_definitions()
def bool_eval(x):
	if not isbool(x):
		print('\n' + logmgs.warning_message('A boolean evaluation was requested for value \'{0}\' but it is not in the bool definitions.').format(x))
		return None
	if x in false_definitions(): return False
	else: return True

# list_intersection: computers intersection between two lists. In case the intersection is null, returns []. {{{2
def list_intersection(list1, list2): return [x for x in list1 if x in list2]

# variable_check is the checking function for variables. It checks for type and range. {{{2
def variable_check(variable, **kwargs):
	if 'custom_message' in kwargs: custom_message = kwargs.get('custom_message')
	else: custom_message = False

	if 'variable_card' in kwargs: variable_card = kwargs.get('variable_card')
	else: variable_card = 'Some variable'
	if 'allowed_type' in kwargs:
		allowed_type = kwargs.get('allowed_type')
		if allowed_type == 'float':
			if not isfloat(variable): raise VariableCheckError(variable_card + ' is not a float.')
		if allowed_type == 'int':
			if not isinteger(variable): raise VariableCheckError(variable_card + ' is not an integer.')
		if allowed_type == 'bool':
			if not isbool(variable): raise VariableCheckError(variable_card + ' is not a boolean.')

	else: allowed_type = 'string'

	if 'allowed_range' in kwargs:
		allowed_range = kwargs.get('allowed_range')
		if allowed_type != 'float':
			print('\n' + logmgs.warning_message('A range check was requested for ' + variable_card + ' but it is not a float, hence the check will not be performed.'), end = '')
		else:
			if allowed_range == '++' and float(variable) <= 0 : raise VariableCheckError(custom_message if custom_message else variable_card + ' has declared value \'{0}\' but should be a positive.'.format(variable))
			if allowed_range == '+0' and float(variable) < 0: raise VariableCheckError(custom_message if custom_message else variable_card + ' has declared value \'{0}\' but should be a non-negative.'.format(variable))
			if allowed_range == '--' and float(variable) >= 0 : raise VariableCheckError(custom_message if custom_message else variable_card + ' has declared value \'{0}\' but should be a negative.'.format(variable))
			if allowed_range == '-0' and float(variable) > 0: raise VariableCheckError(custom_message if custom_message else variable_card + ' has declared value \'{0}\' but should be a non-positive.'.format(variable))
			if allowed_range == '01' and (float(variable) > 1 or float(variable) < 0): raise VariableCheckError(custom_message if custom_message else variable_card + ' has declared value \'{0}\' but should be a positive between 0 and 1.'.format(variable))

	if 'special_floats' in kwargs: special_floats = kwargs.get('special_floats')
	else: special_floats = []

	# The idea behind the list intersection is that if the special_floats keyword argument is passed with any of the definitions in nan_defs(), posinf_defs() or neginf_defs() then the program will compare the variable against all the definitions as well;
	#	for instance, if special_floats = inf, then the function will also compare it against '+inf', 'posinf', 'INF' etc.
	if (	( list_intersection(special_floats, nan_definitions()) == [] and variable in nan_definitions() ) or
		( list_intersection(special_floats, posinf_definitions()) == [] and variable in posinf_definitions() ) or
		( list_intersection(special_floats, neginf_definitions()) == [] and variable in neginf_definitions() )	):
		raise VariableCheckError(custom_message if custom_message else variable_card + ' defined as \'{0}\'.'.format(variable))


	# Used for string or char-type variables check
	if 'allowed_values' in kwargs:
		allowed_values = kwargs.get('allowed_values')
		if variable not in allowed_values: raise VariableCheckError(custom_message if custom_message else variable_card + ' declared value \'{0}\' is not in the allowed roster of values {1}.'.format(variable, allowed_values))

	return True

# data_card_search is a function that takes two arguments: a string 'data_card' and a file object 'file_data'. The function searches the file for the data_card string, sweeping the lines and separating them according to the \t character.  {{{2
# The while loop breaks when the first block of the line is the dataCard or the data card is not found, in which case a DataCardError exception is raised and execution is halted.
# If the data card was found, the function returns the pointer to the line where the dataCard was found.
def data_card_search(data_card,file):
	# Searching for data card
	file.seek(0)
	line = file.readline().strip().split('\t')
	while data_card not in line[0]:
		if line[0] == '':
			print('\n' + logmgs.error_message('data card \'{0}\' was not found!'.format(data_card)))
			exit()

		line = file.readline().strip().split('\t')

	return line

# PARSING SYSTEM DATA ------------------------------------------------------------------------- {{{1
def parse_system_data(file_name):
	with open(file_name,'r') as file:

		print(logmgs.good_message('Loading system data from file \'{0}\'... ').format(file.name), end = '')

		# Serching for generator data finish data card
		data_card_search('END OF SYSTEM DATA', file)
		# Searching for generator data start data card
		line = data_card_search('SYSTEM DATA FOLLOWS',file)

		line = file.readline().strip().split('\t')
		if line[0] != 'CASE ID:': raise DataCardError('CASE ID')
		case_id = line[1]

		line = file.readline().strip().split('\t')
		if line[0] != 'MVA BASE:': raise DataCardError('MVA BASE')
		variable_check(line[1], allowed_type = 'float', allowed_range = '++', variable_card = 'Case MVA BASE base value')
		Sb = float(line[1])

		line = file.readline().strip().split('\t')
		if line[0] != 'VOLTAGE BASE:': raise DataCardError('VOLTAGE BASE')
		if not isfloat(line[1]): raise VariableCheckError('Case voltage base value must be float')
		Vb = float(line[1])

		line = file.readline().strip().split('\t')
		if line[0] != 'DROOP TYPE:': raise DataCardError('DROOP TYPE')
		variable_check(line[1], variable_card = 'Case DROOP TYPE', allowed_values = ['NONE','PERFECT','LINEAR','HYPERBOLIC'])#: raise VariableCheckError(' --> Case Droop is defined as \'{0}\' but must be \'NONE\', \'LINEAR\', \'PERFECT\' or \'HYPERBOLIC\''.format(line[1]))
		droop_type = line[1]

		print(logmgs.green_string(' Done.'))

	return [case_id, Sb, Vb, droop_type]

# PARSING BUS DATA ---------------------------------------------------------------------------- {{{1
def parse_bus_data(file_name):
	with open(file_name,'r') as file:
		print(logmgs.good_message('Loading bus data from file \'{0}\'... ').format(file.name), end = '')
		# Serching for generator data finish data card
		data_card_search('END OF BUS DATA', file)
		# Searching for generator data start data card
		line = data_card_search('BUS DATA FOLLOWS',file)
		bus_list = []
		while True:
			line = file.readline().strip().split('\t')
			line = [x for x in line if x] # Removing empty ''  elements (double tabs)
			
			# Searching for hifen separator: if the line contains 10 hifens, consider bus data over
			if 'END OF BUS DATA' in line[0] : break

			if len(line) != 10:
				print('\n' + logmgs.error_message('bus data requires 10 parameters but {0} are given in row {1}.'.format(len(line), len(bus_list) + 1)))
				exit()

			else:
				# CHECKING BUS NAME ------------------------------------------------
				variable_check(line[0], custom_message = 'There is a bus with name \'NaN\' at position {0}.'.format(len(bus_list) + 1))
				bus_name = line[0]

				# Checking if there already is not another bus with the same name
				for bus in bus_list:
					if bus.name == line[0]: raise VariableCheckError('There are two buses with the same name \'{0}\'!'.format(line[0]))
		
				# CHECKING BUS TYPE ------------------------------------------------
				variable_check(line[1], allowed_values = ['PQ', 'VT', 'PV', 'PdQ', 'dPQ'], variable_card = 'bus \'{0}\' bus type'.format(line[0]), custom_message = 'Bus type of bus \'{0}\' not recognized!'.format(line[0]))
				bus_type = line[1]

				# CHECKING FINAL VOLTAGE -------------------------------------------
				variable_check(line[2], allowed_type = 'float', allowed_range = '++', variable_card = 'bus \'{0}\' voltage magnitude'.format(line[0]))
				final_voltage = float(line[2])

				# CHECKING FINAL ANGLE ---------------------------------------------
				variable_check(line[3], allowed_type = 'float', variable_card = 'bus \'{0}\' voltage angle'.format(line[0]))
				final_angle = float(line[3])

				# CHECKING ACTIVE LOAD ---------------------------------------------
				variable_check(line[4], allowed_type = 'float', allowed_range = '+0', variable_card = 'active load on bus \'{0}\''.format(line[0]))
				pLoad = float(line[4])

				# CHECKING REACTIVE LOAD -------------------------------------------
				variable_check(line[5], allowed_type = 'float', variable_card = 'reactive load on bus \'{0}\''.format(line[0]))
				qLoad = float(line[5])

				# CHECKING SHUNT CONDUCTANCE gsh -----------------------------------
				variable_check(line[6], allowed_type = 'float', allowed_range = '+0', variable_card = 'reactive load on bus \'{0}\''.format(line[0]))
				gsh = float(line[6])

				# CHECKING SHUNT SUSCEPTANCE bsh -----------------------------------
				variable_check(line[7], allowed_type = 'float', allowed_range = '+0', variable_card = 'shunt susceptance on bus \'{0}\''.format(line[0]))
				bsh = float(line[7])

				# CHECKING BUS VOLTAGE LIMITS --------------------------------------
				variable_check(line[8], allowed_type = 'float', allowed_range = '++', special_floats = ['inf'], variable_card = 'maximum voltage magnitude of bus \'{0}\''.format(line[0]))
				vmax = float(line[8])

				variable_check(line[9], allowed_type = 'float', allowed_range = '+0', variable_card = 'minimum voltage magnitude of bus \'{0}\''.format(line[0]))
				vmin = float(line[9])

				if vmax <= vmin:
					print(logmgs.error_message('Minimum voltage limit of bus \'{0}\' is greater or equal than its maximum voltage limit.'.format(line[0])))
					exit()

				# CREATING BUS INSTANCE --------------------------------------------
				bus_list.append(cL.bus( 0, bus_name, bus_type, final_voltage, final_angle, pLoad, qLoad, gsh, bsh, vmax, vmin))	

		# Checking if there is a VT bus
		for bus in bus_list:
			if bus.bus_type == 'VT': break
		else:
			print('\n' + logmgs.error_message('no VT bus was found. At least one must be present in the system.'))
			exit()

		# i is the bus number counter. It is used to assign the bus numbers that will be used by the program; bus numbers are assigned in the order they appear in the netfile.
		for i in range(len(bus_list)): bus_list[i].number = i

		print(logmgs.green_string(' Done.'))

	counter = 0
	for bus in bus_list:
		bus.number = counter
		counter += 1

	return bus_list

# PARSING BRANCH DATA ------------------------------------------------------------------------- {{{1
def parse_branch_data(file_name, bus_list):
	with open(file_name,'r') as file:
		print(logmgs.good_message('Loading branch data from file \'{0}\'... ').format(file.name), end = '')

		# Serching for branch data finish data card
		data_card_search('END OF BRANCH DATA', file)
		# Searching for branch data start data card
		line = data_card_search('BRANCH DATA FOLLOWS',file)

		branch_list = []
		while True:
			line = file.readline().strip().split('\t')
			line = [x for x in line if x] # Removing empty ''  elements (double tabs)
			
			# Searching for hifen separator: if the line contains 10 hifens, consider bus data over
			if 'END OF BRANCH DATA' in line[0] : break
			if len(line) != 12:
				print('\n' + logmgs.error_message('branch data requires 12 parameters but {0} are given in row {1}.'.format(len(line), len(branch_list) + 1)))
				exit()

			else:
				# Checking tapbus and loadbus -------------------------------------
				from_bus = line[0]
				to_bus = line[1]
				# Testing if from_bus and to_bus are the same
				if from_bus == to_bus: 
					print('\n' + logmgs.error_message('there is a branch that has coinciding tap bus and load bus (\'{0}\')'.format(from_bus)))
					exit()
				# Testing if to_bus and from_bus are defined in the bus list
				from_bus_flag = False
				to_bus_flag = False
				for bus in bus_list:
					if bus.name == from_bus: from_bus_flag = True
					if bus.name == to_bus: to_bus_flag = True

				if not from_bus_flag:
					print('\n' + logmgs.error_message('branch with tap bus \'{0}\' is invalid because such bus is not defined in the bus list.'.format(from_bus)))
					exit()
				if not to_bus_flag:
					print('\n' + logmgs.error_message('branch with load bus \'{0}\' is invalid because such bus is not defined in the bus list.'.format(to_bus)))
					exit()


				# CHECKING BRANCH TYPE ---------------------------------------------
				variable_check(line[2], allowed_values = ['F', 'VTV', 'VTQ', 'VP'], variable_card = 'branch from bus \'{0}\' to bus \'{1}\''.format(from_bus, to_bus))
				branch_type = line[2]
				
				# CHECKING R AND X -------------------------------------------------
				variable_check(line[3], allowed_type = 'float', allowed_range = '+0', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' resistance'.format(from_bus, to_bus))
				variable_check(line[4], allowed_type = 'float', allowed_range = '++', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' reactance'.format(from_bus, to_bus))
				r = float(line[3])
				x = float(line[4])

				# CHECKING GSH AND BSH ---------------------------------------------
				variable_check(line[5], allowed_type = 'float', allowed_range = '+0', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' shunt conductance'.format(from_bus, to_bus))
				variable_check(line[6], allowed_type = 'float', allowed_range = '+0', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' shunt susceptance'.format(from_bus, to_bus))
				gsh = float(line[5])
				bsh = float(line[6])

				# CHECKING TRANSFORMER RATION, PHASE AND MAX/MIN VALUES ------------
				variable_check(line[7], allowed_type = 'float', allowed_range = '01', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer ratio'.format(from_bus, to_bus))
				variable_check(line[8], allowed_type = 'float', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer phase angle'.format(from_bus, to_bus))
				a = float(line[7])
				phi = float(line[8])*np.pi/180

				if branch_type in ['VTV', 'VTQ']:
					variable_check(line[9], allowed_type = 'float', allowed_range = '01', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer minimum tap'.format(from_bus, to_bus))
					variable_check(line[10], allowed_type = 'float', allowed_range = '01', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer maximum tap'.format(from_bus, to_bus))
					mintp = float(line[9]) 
					maxtp = float(line[10])
				elif branch_type == 'VP':
					variable_check(line[9], allowed_type = 'float', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer minimum phase angle'.format(from_bus, to_bus))
					variable_check(line[10], allowed_type = 'float', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer maximum phase angle'.format(from_bus, to_bus))
					mintp = float(line[9])*np.pi/180
					maxtp = float(line[10])*np.pi/180
				elif branch_type == 'F':
					variable_check(line[9], allowed_type = 'float', special_floats = ['nan'], variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer minimum tap/phase'.format(from_bus, to_bus))
					variable_check(line[10], allowed_type = 'float', special_floats = ['nan'], variable_card = 'branch from bus \'{0}\' to bus \'{1}\' transformer maximum tap/phase'.format(from_bus, to_bus))
					mintp = float(line[9])
					maxtp = float(line[10])
				if mintp >= maxtp:
					print('\n' + logmgs.error_message('branch from bus \'{0}\' to bus \'{1}\' cannot have a maximum transformer tap/phase less than or equal than the minimum.'.format(from_bus, to_bus)))
					exit()

				# CHECKING BRANCH STATUS -------------------------------------------
				variable_check(line[11], allowed_type = 'bool', variable_card = 'branch from bus \'{0}\' to bus \'{1}\' status'.format(from_bus, to_bus))
				status = bool_eval(line[11])
				new_branch = cL.branch( 0, from_bus, to_bus, r, x, gsh, bsh, a, phi, mintp, maxtp, status)
				branch_list.append(new_branch)

		print(logmgs.green_string(' Done.'))
		counter = 0
		for branch in branch_list:
			branch.number = counter
			counter += 1

		return branch_list

# PARSING GENERATOR DATA ---------------------------------------------------------------------- {{{1
def parse_generator_data(gendata_file, bus_list, system_data):
	with open(gendata_file,'r') as file:
		print(logmgs.good_message('Loading generator data from file \'{0}\'... ').format(file.name), end = '')
		# Serching for generator data finish data card
		data_card_search('END OF GENERATOR DATA', file)
		# Searching for generator data start data card
		line = data_card_search('GENERATOR DATA FOLLOWS',file)

		# Extracting system droop type from system data
		droop_type = system_data[3]
		gen_list = []
		while True:
			line = file.readline().strip().split('\t')
			line = [x for x in line if x] # Removing empty ''  elements (double tabs)
			
			# Searching for hifen separator: if the line contains 10 hifens, consider bus data over
			if 'END OF GENERATOR DATA' in line[0] : break
			# Checking if parameters are missing or surplus
			if len(line) != 13:
				print('\n' + logmgs.error_message('generator data requires 13 parameters but {0} are given in row {1}.'.format(len(line), len(gen_list) + 1)))
				exit()
			else:
				# Checking if the generator is attached to a bus
				for bus in bus_list:
					if str(line[0]) == bus.name:
						break
				else:	# If break was not done it was because the tap bus name declared does not exist
					print('\n' + logmgs.error_message('there is a generator attached to a non-declared bus \'{0}\'.'.format(line[0])))
					exit()
			
				# Checking if the generator is not attached to a bus that already has a generator
				for gen in gen_list:
					if str(line[0]) == gen.bus_name:
						print('\n' + logmgs.error_message('there are two generators attached to bus \'{0}\''.format(line[0])))
						exit()

				# Finally, go through customary checks
				bus_name = line[0] 

				# Checking pgen and qgen ------------------------------------------
				variable_check(line[1], allowed_type = 'float', allowed_range = '+0', variable_card = 'active power generation of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[2], allowed_type = 'float', variable_card = 'reactive power generation of generator at bus \'{0}\''.format(bus_name))
				p_gen = float(line[1])
				q_gen = float(line[2])

				# Checking P and Q limits ------------------------------------------
				variable_check(line[3], allowed_type = 'float', special_floats = ['inf'], variable_card = 'maximum reactive power generation of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[4], allowed_type = 'float', special_floats = ['-inf'], variable_card = 'minimum reactive power generation of generator at bus \'{0}\''.format(bus_name))
				q_max = float(line[3])
				q_min = float(line[4])

				variable_check(line[5], allowed_type = 'float', allowed_range = '+0', special_floats = ['inf'], variable_card = 'maximum active power generation of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[6], allowed_type = 'float', allowed_range = '+0', variable_card = 'minimum active power generation of generator at bus \'{0}\''.format(bus_name))
				p_max = float(line[3])
				p_min = float(line[4])


				# Checking rated values --------------------------------------------
				variable_check(line[7], allowed_type = 'float', allowed_range = '++', variable_card = 'voltage reference of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[8], allowed_type = 'float', allowed_range = '++', variable_card = 'rated active power of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[9], allowed_type = 'float', allowed_range = '++', variable_card = 'rated reactive power of generator at bus \'{0}\''.format(bus_name))
				variable_check(line[10], allowed_type = 'float', allowed_range = '++', variable_card = 'rated voltage of generator at bus \'{0}\''.format(bus_name))
				v_ref = float(line[7])
				p_rated = float(line[8])
				q_rated = float(line[9])
				v_rated = float(line[10])

				# Checking Droop ramp values ---------------------------------------
				variable_check(line[11], allowed_type = 'float', allowed_range = '--', special_floats = ['nan'], variable_card = 'active Droop ramp of generator at bus \'{0}\''.format(bus_name))
				p_ramp = float(line[11])
				variable_check(line[12], allowed_type = 'float', allowed_range = '--', special_floats = ['nan'], variable_card = 'reactive Droop ramp of generator at bus \'{0}\''.format(bus_name))
				q_ramp = float(line[12])

				if droop_type != 'NONE' and cL.get_bus(bus_name, bus_list).bus_type == 'dPQ' and np.isnan(p_ramp):
					print('\n' + logmgs.error_message('generator in bus \'{0}\' has a \'NaN\' active Droop ramp yet the case Droop type is not \'NONE\' and the generator bus is defined as \'dPQ\' type. Either enter a valid (negative) value for the ramp or define the bus as PQ or PV.'.format(line[0])))
					exit()

				if droop_type != 'NONE' \
				and ( cL.get_bus(bus_name, bus_list).bus_type == 'dPQ' or cL.get_bus(bus_name, bus_list).bus_type == 'PdQ' ) \
				and np.isnan(p_ramp) :
					print('\n' + logmgs.error_message('generator in bus \'{0}\' has a \'NaN\' reactive Droop ramp yet the case Droop type is not \'NONE\' and the generator bus is defined as \'{1}\' type. Either enter a valid (negative) value for the ramp or define the bus as PQ.'.format(bus_name, cL.get_bus(bus_name, bus_list).bus_type )))
					exit()
				

				# Creating the new generator instance to be added
				new_gen = cL.generator(bus_name, p_gen, q_gen, q_max, q_min ,p_max, p_min, v_ref, p_rated, q_rated, v_rated, p_ramp, q_ramp)
				# If the parameters of the generator were given in relation with the generator PU system, they should be converted to the system's PU			
				gen_list.append(new_gen)

			#	print(newGen)

		# Checking if all Droop-type buses are attached to a generator -----
		for bus in bus_list:
			if bus.bus_type in ['dPQ', 'PdQ']:
				if not cL.is_gen(bus.name, gen_list): raise DataCardError(' --> Bus \'{0}\' is defined as {1} but does not have a generator attached.'.format(bus.name, bus.bus_type))				

				
		print(logmgs.green_string(' Done.'))
		return gen_list

# PARSING MACHINE DATA ------------------------------------------------------------------------ {{{1
def parse_machine_data(machinedata_file, gen_list): return []

# PARSING FAULT DATA -------------------------------------------------------------------------- {{{1
def parse_fault_data(faultdata_file, branch_list, bus_list): return []
	
# MAIN FUNCTION ------------------------------------------------------------------------------- {{{1
def load_case(net_file,**kwargs):
	if 'systemdata_file' in kwargs: systemdata_file = kwargs.get('systemdata_file')
	else: systemdata_file = net_file

	if 'busdata_file' in kwargs: busdata_file = kwargs.get('busdata_file')
	else: busdata_file = net_file

	if 'branchdata_file' in kwargs: branchdata_file = kwargs.get('branchdata_file')
	else: branchdata_file = net_file

	if 'gendata_file' in kwargs: gendata_file = kwargs.get('gendata_file')
	else: gendata_file = net_file

	if 'machinedata_file' in kwargs: machinedata_file = kwargs.get('machinedata_file')
	else: machinedata_file = net_file

	if 'faultdata_file' in kwargs: faultdata_file = kwargs.get('faultdata_file')
	else: faultdata_file = net_file

	print(logmgs.header_message('LOADING NET FILE \'{0}\''.format(net_file)))
	system_data = parse_system_data(systemdata_file)
	bus_list = parse_bus_data(busdata_file)
	branch_list = parse_branch_data(branchdata_file, bus_list)
	gen_list = parse_generator_data(gendata_file, bus_list, system_data)
	machine_list = parse_machine_data(machinedata_file, gen_list)
	fault_list = parse_fault_data(faultdata_file, branch_list, bus_list)
	
	case = cL.case(system_data, bus_list, branch_list, gen_list, machine_list, fault_list)
	print(logmgs.green_string(' -> Case\'{0}\' from file \'{1}\' loaded successfully. No errors reported. \n'.format(system_data[0], net_file)))
	#case.update_matrixes()

	return case
