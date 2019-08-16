# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: arguments library for pyEPSim
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.0
# DESCRIPTION: this file describes the arguments that can be passed to pyEPSim.
# -------------------------------------------------

from libraries import misc_functions as miscF
clear = miscF.clear

import argparse
parser = argparse.ArgumentParser(
	prog = 'pyEPSim Python Electric Power System Simulator',
	description = 'Dynamical simulator that takes a system data in the IEEE CDF and a measurements file in the LACO CDF.',
	)

def parseArguments():
	# Arguments and options
	#	Creates arguments:
	#	--net [NETWORK FILE] is the target *.net system network file. This string is stored in a 'networkFile' variable. The network file describes the system topologically and parametrically: how many buses, which buses are connected to which, how many generators and how they are connected; the parameters of these generators, the faults applicable and their quality.
	parser.add_argument(	"net",
				type = str,
				metavar = 'NETFILE',
				help='Target network file')

	# VERBOSE option
	parser.add_argument(	'--verbose', '-v',
				type = int,
				choices = [0,1,2],
				default = [0],
				nargs = 1,
				help = 'Verbose level. 0 for no output, 1 for critical output, 2 for detailed output. Default is 0.')

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

	return args
