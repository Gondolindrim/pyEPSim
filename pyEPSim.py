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

# Numpy printopt so that the console does not limit the output column number
import numpy as np
np.set_printoptions(linewidth=999)

# Loading libraries
import libraries.classes as cL		# classes.py contains the bus, branch, generator and case classes
import libraries.arguments as pA	# arguments.py contains the available arguments passed to the program
import libraries.loadCase as lC		# loadCase.py contains the case-loading routine which reads the net file
import libraries.powerFlow as pF	# powerFlow.py contains the Newton-Raphson method used for power flow calculation.

# Importing tabulate dependency for pretty tabular prints
from libraries.tabulate.tabulate import tabulate

# -------------------------------------------------

# -------------------------------------------------
# (2) DESCRIBRING SYSTEM: READING FILE AND LOADING CASE {{{1
# ----------------------------------------------------------

args = pA.parseArguments()
verbose = args.verbose
netFile = args.net

case = lC.loadCase(netFile)
case.runPowerFlow()

print(case)
