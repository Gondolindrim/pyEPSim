# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: miscellaneous functions library for pyEPSim
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.0
# DESCRIPTION: this library contains miscellaneous functions that don't fit into numerical or computational categories.
# -------------------------------------------------

# (0) OS and SYS libraries for defining clear terminal command
import os
import sys

import numpy as np

# (1) Set printout options: this allows the program printout to span the whole terminal
def setPrintoptions(): np.set_printoptions(suppress=False,linewidth=999,threshold=999)

# (2) Clear  function: clears the terminal. Similar to MATLAB's cls and bash's cls
def clear(): os.system('cls' if os.name == 'nt' else 'clear')

# (3) "pause" function: halts program execution until ENTER key is typed
def pause(string): programPause = input(string)

# (4) readtabline function
# readtabline reads a line in a given file f then arranges it in an array which vectors are the strings
#	in the line that are separated by the tab character '\t'. It also removes the new line '\n' character
#	from the last slot in the line.
def readtabline(f): return f.readline().strip().split('\t')
