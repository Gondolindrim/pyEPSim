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
from libraries.dynamicSimulation import dynamicSimulation as dS
# Importing tabulate dependency for pretty tabular prints
from libraries.tabulate.tabulate import tabulate

from copy import deepcopy as copy
# -------------------------------------------------

# -------------------------------------------------
# (2) DESCRIBRING SYSTEM: READING FILE AND LOADING CASE {{{1
# ----------------------------------------------------------
np.set_printoptions(threshold=np.inf)
args = pA.parseArguments()
verbose = args.verbose
net_file = args.net

case = lC.load_case(net_file)
success = case.run_power_flow()
if success: print(case)
#print(case)
#Yred, rCase = case.reduceMatrixes()
#
## TESTING REDUCTION ALGORITHM
#
## TEST (1): rCase and case power flow should be the exact same
#rCase.runPowerFlow()
#print(case)
#print(rCase)
#
## TEST (2): currents in case should be the same, be it calculated through power flux or conductance matrix
#V = [case.busData[k].finalVoltage for k in range(case.nBus)]
#theta = [np.pi*case.busData[k].finalAngle/180 for k in range(case.nBus)]
#
#P = np.array([bus.pGen for bus in case.busData])/case.Sb
#Q = np.array([bus.qGen for bus in case.busData])/case.Sb
#V = [V[k]*np.e**(1j*theta[k]) for k in range(case.nBus)]
#YL = np.diag([ (case.busData[k].pLoad  - 1j*case.busData[k].qLoad)/case.busData[k].finalVoltage**2 for k in range(case.nBus)])/case.Sb
#Y = copy(case.Y) + YL
#IC = [ np.conj((P[k] + 1j*Q[k])/(V[k])) for k in range(case.nBus)]	# Currents calculated through power flow
#I = Y @ V	# Currents calculated throgh conductance matrix
#print(' ---> TEST 2: calculating norm(IC - I) = {}'.format(np.linalg.norm(IC - I)))
#
## TEST (3): currents in rcase should be the same, be it calculated through power flux or conductance matrix
#P = np.array([bus.pGen for bus in rCase.busData])/rCase.Sb
#Q = np.array([bus.qGen for bus in rCase.busData])/rCase.Sb
#V = [rCase.busData[k].finalVoltage for k in range(rCase.nBus)]
#theta = [np.pi*rCase.busData[k].finalAngle/180 for k in range(rCase.nBus)]
#V = [V[k]*np.e**(1j*theta[k]) for k in range(rCase.nBus)]
#
#YL = np.diag([ (rCase.busData[k].pLoad  - 1j*rCase.busData[k].qLoad)/rCase.busData[k].finalVoltage**2 for k in range(rCase.nBus)])/rCase.Sb
#Y = copy(rCase.Y) + YL
#IC = [ np.conj((P[k] + 1j*Q[k])/(V[k])) for k in range(rCase.nBus)]
#I = Y @ V
#print(' ---> TEST 3: calculating norm(IC - I) = {}'.format(np.linalg.norm(IC - I)))
#
## TEST (4): currents in rCase generators should be the same when calculated through power flux or Yred matrix
#P = np.array([bus.pGen for bus in rCase.busData])/rCase.Sb
#Q = np.array([bus.qGen for bus in rCase.busData])/rCase.Sb
#V = [rCase.busData[k].finalVoltage for k in range(rCase.nBus)]
#theta = [np.pi*rCase.busData[k].finalAngle/180 for k in range(rCase.nBus)]
#V = [V[k]*np.e**(1j*theta[k]) for k in range(rCase.nBus)]
#
#YL = np.diag([ (rCase.busData[k].pLoad  - 1j*rCase.busData[k].qLoad)/rCase.busData[k].finalVoltage**2 for k in range(rCase.nBus)])/rCase.Sb
#Y = copy(rCase.Y) + YL
#IC = [ np.conj((P[k] + 1j*Q[k])/(V[k])) for k in range(rCase.nGen)]
#I = Yred @ [V[k] for k in range(rCase.nGen)]
#print(' ---> TEST 4: calculating norm(IC - I) = {}'.format(np.linalg.norm(IC - I)))
#print(' ---> I = {}'.format(I))
#
# TEST (5): dynamic simulation
#disturbanceData = ['Bus003', 2*(10 + 1j*5), 1]
#dS(case, disturbanceData, 40)
