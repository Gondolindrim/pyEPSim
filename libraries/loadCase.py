# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: loadCase function
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.2
# DESCRIPTION: # loadcase is a function that takes the pointer "fileName" to a file and returns a case structure (see (2.6)).
# The structure is composed of bus, branch, generator and fault data, comprised of the lists of
# buses, branches, generators and faults. Each of these have a particular set of attributes that
# are explained in their respective object definitions.
# -------------------------------------------------

from libraries import classes as cL

class dataCardError(Exception):
	def __init__(self, msg):
		self.code = msg

def loadCase(fileName):
	print(' --> Loading case file \'{0}\'...'.format(fileName))
	with open(fileName,'r') as fileData:

		line = [1]
		while True:
			if line[0] == 'MVA base:':
				line = fileData.readline().strip().split('\t')
				Sb = float(line[1])	# Extracting power base value
				break
			elif line == '':
				break
				raise dataCardError(' --> Netfile error: \'MVA base\' data card was not found! Please define a base power value in the netfile.')

		while True:
			if line[0] == 'Voltage base:':
				line = fileData.readline().strip().split('\t')
				Vb = float(line[1])	# Extracting voltage base value
			elEOFError(' --> Netfile error: \'Voltage base\' card was not found! Please define a base voltage value in the netfile.')

		# Bus data parameters ------------------------------------------
		# Searching for bus data start card
		while line[0] != 'BUS DATA FOLLOWS': line = fileData.readline().strip().split('\t')

		# Acquiring bus data
		busList = []
		while True:
			line = fileData.readline().strip().split('\t')
			if line[0] == '-999': break
			else: busList.append(cL.bus( line[0], line[1], line[7], line[8], line[9], line[10], line[15], line[16]))
		else: raise EOFError(' --> Netfile error: the endcard '-999' was not found when acquiring bus data')
			
		# Branch data parameters ---------------------------------------
		# Searching for branch data start card
		while line[0] != 'BRANCH DATA FOLLOWS':	line = fileData.readline().strip().split('\t')

		# Acquiring branch data
		branchList = []	
		while line[0] != '-999': branchList.append(cL.branch(line[0],line[1],line[6],line[7],line[8],1/float(line[14])))

		# Generator data parameters ------------------------------------
		# Searching for generator data start card
		while line[0]!= 'GENERATOR DATA FOLLOWS': line = fileData.readline().strip().split('\t')
		fileData.readline()	# Skip the '---' line	

		genList = []	# genList was named like so because there already is a genData variable in the program
		while line[0] != '-999':
			line = fileData.readline().strip().split('\t')
			genList.append(cL.generator(line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16]))
		
		# Fault data ---------------------------------------------------
		# Searching for fault data start card 
		while line[0] != 'FAULT DATA FOLLOWS': line = fileData.readline().strip().split('\t')
		next(fileData)	# Skip the '---' line	

		faultList = []	# As the same case with genList, faultList was named like so because there already is a faultData variable in the program
		while line[0] != '-999':
			line = fileData.readline().strip().split('\t')
			faultList.append(cL.fault( line[0], line[1], line[2]))

		print(' Done.')

	return case(busList,branchList,genList,faultList,Sb,Vb)
