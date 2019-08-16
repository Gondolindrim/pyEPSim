#!/usr/bin/python
# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: pyEPSim, a Python dynamical simulator for Electric Power Systems.
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: V5.1.1
# DESCRIPTION: this program is a power flow calculator for an electric power system. It receives
# BIBLIOGRAPHY:
# [1] 
# [2]
# [3]
# -------------------------------------------------

__author__ = "Alvaro Augusto Volpato and Luís Fernando Costa Alberto"
__license__ = "GPL"
__version__ = "5.2"
__maintainer__ = "Alvaro Augusto Volpato"
__email__ = "alvaro.volpato@usp.br"
__status__ = "Development"

# -------------------------------------------------
# (1) IMPORTING LIBRARIES  {{{1
# -------------------------------------------------

import numpy as np
desired_width = 999
np.set_printoptions(linewidth=desired_width)

# (1.8) Loading 
import libraries.classes as cL		# classes.py contains the bus, branch, generator and case classes
import libraries.arguments as pA	# arguments.py contains the available arguments passed to the program
import libraries.loadCase as lC		# loadCase.py contains the case-loading routine which reads the net file
import libraries.powerFlow as pF	# powerFlow.py contains the Newton-Raphson method used for power flow calculation.

from libraries.tabulate.tabulate import tabulate

from libraries.color import color as color

# -------------------------------------------------

# -------------------------------------------------
# (2) DESCRIBRING SYSTEM: READING FILE AND LOADING CASE {{{1
# ----------------------------------------------------------

args = pA.parseArguments()
verbose = args.verbose
netFile = args.net

case14 = lC.loadCase(netFile)

print(case14)
	
print(' --> Beggining power flow...')
V, theta, r, elapsed, itCount, success = pF.powerFlow(case14,deltaMax=100)
print(' Done.\n')

tabHeaders = ['Name','Number','Voltage (p.u.)','Angle (rad)']

tabRows = []
for bus in case14.busData:
	tabRows.append([bus.name, bus.number, V[bus.number-1], theta[bus.number-1]])
	
resultsTable = tabulate(tabRows,headers=tabHeaders)

# (6.13) Printing convergence results
if success:
	print(' --> Power flow result:')
	print(resultsTable)
#	if verbose > 1: print('\n --> Residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,r))

else: print(' --> Power flow method not successful.')
