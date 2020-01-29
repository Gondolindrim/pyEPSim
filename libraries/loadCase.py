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

import numpy as np

import libraries.classes as cL

from copy import deepcopy
copy = deepcopy

# dataCardError is a custom-created exception raised whenever a start or end card is not found.
class dataCardError(Exception):
	pass

# dataCardSearch is a function that takes two arguments: a string 'dataCard' and a file object 'fileData'. The function searches the file for the dataCard string, sweeping the lines and separating them according to the \t character. 
# The while loop breaks when the first block of the line is the dataCard or the data card is not found, in which case a dataCardError exception is raised and execution is halted.
# If the data card was found, the function returns the pointer to the line where the dataCard was found.
def dataCardSearch(dataCard,fileData):
	while True:
		line = fileData.readline().strip().split('\t')
		if line[0] == dataCard: break
		elif line[0] == '':
			raise dataCardError(' --> Netfile error: \'{0}\' card was not found!'.format(dataCard))
			break
	return line

def loadCase(fileName,**kwargs):
	print(' --> Loading case file \'{0}\'...'.format(fileName),end='')

	caseID = 'Case {0}'.format(fileName)

	with open(fileName,'r') as fileData:

		line = [1]
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == 'Case id:':
				caseID = line[1]	# Extracting power base value
				break
			elif line[0] == '':
				break

	with open(fileName,'r') as fileData:

		line = [1]
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == 'MVA base:':
				Sb = float(line[1])	# Extracting power base value
				break
			elif line[0] == '':
				raise dataCardError(' --> Netfile error: \'MVA base\' data card was not found! Please define a base power value in the netfile.')
				break

		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == 'Voltage base:':
				Vb = float(line[1])	# Extracting voltage base value
				break
			elif line[0] == '':
				raise dataCardError(' --> Netfile error: \'Voltage base\' card was not found! Please define a base voltage value in the netfile.')
				break

		# Bus data parameters ------------------------------------------
		# Searching for bus data start card
		line = dataCardSearch('BUS DATA FOLLOWS',fileData)
		line = fileData.readline()	# Skip --- line

		# Acquiring bus data
		busList = []
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == '-999': break
			elif line[0] == '':
				raise dataCardError(' >> Netfile error: the endcard \'-999\' was not found when acquiring bus data')
				break
			else:
				busList.append(cL.bus( 0, line[0], line[3], line[6], line[7], line[8], line[9], line[10], line[15],line[4],line[5]))
		
		# Reorganizing buses so the first bus is the slack
		for i in range(len(busList)):
			if busList[i].PVtype == 'VT':
				busList[i], busList[0] = busList[0], busList[i]
				break
		else:
			raise dataCardError (' >> Netfile error: no VT bus was found! Please assign a reference bus.')

		# i is the bus number counter. It is used to assign the bus numbers that will be used by the program; bus numbers are assigned in the order they appear in the netfile.
		for i in range(len(busList)): busList[i].number = i
			
		# Branch data parameters ---------------------------------------
		# Searching for branch data start card
		line = dataCardSearch('BRANCH DATA FOLLOWS',fileData)
		line = fileData.readline()	# Skip --- line

		branchList = []	
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == '-999': break
			elif line[0] == '':
				raise dataCardError(' >> Netfile error: the endcard \'-999\' was not found when acquiring branch data')
				break
			else:
				# Searching for the number of the tap bus
				for bus in busList:
					if str(line[0]) == bus.name:
						fromBus = bus.number
						break
				else:	# If break was not done it was because the tap bus name declared does not exist
					raise dataCardError(' >> Netfile error: there is a branch declared with non-declared tap bus \'{0}\''.format(line[0]))
					break
				
				# Searching for the number of the Z bus
				for bus in busList:
					#print(bus.name)
					if str(line[1]) == bus.name:
						toBus = bus.number
						break
				else:	# If break was not done it was because the Z bus name declared does not exist
					raise dataCardError(' >> Netfile error: there is a branch declared with non-declared Z bus \'{0}\''.format(line[1]))
					break

				if float(line[14]) != 0: branchList.append(cL.branch(0, fromBus,toBus,line[6],line[7],line[8],1/float(line[14])))
				else: branchList.append(cL.branch(0,fromBus,toBus,line[6],line[7],line[8],1))

		# Just likle with the buses, i is the branch number counter. It is used to assign the branch numbers that will be used by the program; branch numbers are assigned in the order they appear in the netfile.
		for i in range(len(branchList)): branchList[i].number = i

		# Generator data parameters ------------------------------------
		# Searching for generator data start card
		line = dataCardSearch('GENERATOR DATA FOLLOWS',fileData)
		line = dataCardSearch('PU BASE:',fileData)
		genDataPUReference = str(line[1])
		if genDataPUReference != 'SYSTEM' and genDataPUReference != 'GENERATOR':
			raise dataCardError(' >> Generator PU reference declared is wrong or was not found. Please inform a valid option SYSTEM or GENERATOR.')

		line = fileData.readline()	# Skip --- line

		genList = []	# genList was named like so because there already is a genData variable in the program
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == '-999': break
			elif line[0] == '':
				raise dataCardError(' --> Netfile error: the endcard \'-999\' was not found when acquiring generator data')
				break
			else:
				for bus in busList:
					if str(line[0]) == bus.name:
						genBus = bus.number
						break
				else:	# If break was not done it was because the tap bus name declared does not exist
					raise dataCardError(' >> Netfile error: there is a generator declared with non-declared bus \'{0}\''.format(line[0]))
					break


				# Creating the new generator instance to be added
				newGen = cL.generator(genBus,line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16])

				# If the parameters of the generator were given in relation with the generator PU system, they should be converted to the system's PU
				if genDataPUReference == 'GENERATOR' or genDataPUReference == 'generator':
					newGen.H *= newGen.ratedPower/Sb
					newGen.ra *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xL *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xd *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xPd *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xPPd *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xq *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xPq *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)
					newGen.xPPq *= (newGen.ratedPower/newGen.ratedVoltage**2)/(Sb/Vb**2)

				genList.append(newGen)
	

		# Sorting generators by their bus number
		for i in range(len(genList)):
			for j in range(len(genList)):
				if genList[i].busNumber < genList[j].busNumber: genList[i], genList[j] = genList[j], genList[i]
		
		# Fault data ---------------------------------------------------
		# Searching for fault data start card 
		line = dataCardSearch('FAULT DATA FOLLOWS',fileData)
		line = fileData.readline()	# Skip --- line

		faultList = []	# As the same case with genList, faultList was named like so because there already is a faultData variable in the program

		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == '-999': break
			elif line[0] == '':
				raise dataCardError(' --> Netfile error: the endcard \'-999\' was not found when acquiring fault data')
				break
			else: faultList.append(cL.fault( line[0], line[1], line[2]))			

		print(' Done.')

	case = cL.case(caseID,busList,branchList,genList,faultList,Sb,Vb)
	case.updateMatrixes()

	return case
