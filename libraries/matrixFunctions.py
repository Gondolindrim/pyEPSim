# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: pyEPSim matrix functions library
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.3
# DESCRIPTION: this file contains the matrix functions used to build the jacobian of the Newton-Raphson method of pyEPSim.
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

# Common functions {{{1

# Lambda function: common to all derivatives
def LambdaP(V,a,b,d): return V*( a*cos(d) + b*sin(d) )
#}}}1

# Functions used for state estimation and power flow {{{1
def dPdT(V,theta,case): #{{{2
	H = np.zeros((int(case.isP.sum()), int(case.isT.sum())))	# Note that H has the size nxm, where n is the amount of active power variables to be estimated and m is the amount of angle variables to be estimated. Also note that this sum takes account not only of the bus power injections but also the transmission line power fluxes. This is where the isP.sum() and isT.sum() come from: since the elements of these matrixes are either zero or one, the sum of their elements give precisely the sizes of matrix H. This pattern also repeats itself to the other matrixes that make the jacobian, so this "trick" is employed in all of them.
	i = 0	# Starting counter
	for j in range(case.nBus):
		for n in range(case.nBus):
			if case.isP[j,n]:
				p = 0
				for k in range(case.nBus):
					if case.isT[k]:
						if j==n: # If measuring power injection
							if k==j: # Calculating dPj/dTj. Note that since Pj = -sum(Pkj) for all k connected to j, from whence comes the possibility og using the same dPj_dTk function with the summation.
								H[i,p] = V[j]*sum([LambdaP(V[m],case.B[j,m],-case.G[j,m],theta[k] - theta[m]) for m in case.K[j]])
							elif (k in case.K[j]): # Calculating dPj/dTk when j and k are connected
								H[i,p] = V[j]*LambdaP(V[k],-case.B[j,k],case.G[j,k],theta[j] - theta[k])
							# When j and k are not connected, then dPj/dTk = 0, but since the default value of H[i,p] is already zero, this case is not analyzed

						else:	# If measuring power flow: calculating dPjn/dTk
							# if k!=j and k!=n, since power transfer is function of only angles j and n, this is zero. But since the default value of H[i,p] is already zero, then this case is not alayzed.
							if k==j: H[i,p] = V[j]*LambdaP(V[n],case.B[j,n],-case.G[j,n],theta[k] - theta[n])
							else: H[i,p] = V[j]*LambdaP(V[n],-case.B[j,n],case.G[j,n],theta[k] - theta[n])
						#print(' --> dPdT[{},{}] = dP({},{})/dT({})'.format(i,k,j,n,k))
						p+=1
				i+=1
	# Deleting the reference bar
	#print(H)
	return H

def dPdV(V,theta,case):#{{{2
	N = H = np.zeros((int(case.isP.sum()), int(case.isV.sum())))
	i = 0
	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isP[j,n]:
				p = 0
				for k in range(case.nBus):
					if case.isV[k]:
						if j==n: # Calculating dPj/dTk
							if k==j: N[i,p] = 2*V[j]*case.G[j,j] + sum([LambdaP(V[m],case.G[j,m],case.B[k,m],theta[j] - theta[m]) for m in case.K[j]])
							elif k in case.K[j]: N[i,p] = LambdaP(V[k],case.G[j,k],case.B[j,k],theta[j] - theta[k])
							# The case where k != j,n is not considered as this case equals to zero and that is already the default value of N.
						else: # Calculating dPjn/dTk
							# Calculating dPjn/dPj
							if k == j: N[i,p] = 2**V[j]*case.Gsh[j,n] + LambdaP(V[n],case.G[j,n],case.B[j,n],theta[j] - theta[n])
							else: N[i,p] = LambdaP(V[j],case.G[j,n],case.B[j,n],theta[j] - theta[n])
							# The case where k != j,n is not considered as this case equals to zero and that is already the default value of N.
						#print(' --> dPdV[{},{}] = dP({},{})/dV({}) = {}'.format(i,p,j,n,k,N[i,p]))
						p+=1
				i+=1
	#print(N)
	return N

def dQdT(V,theta,case): #{{{2
	M = np.zeros((int(case.isQ.sum()), int(case.isT.sum())))
	i = 0
	for j in range(case.nBus):
		for n in range(case.nBus):
			if case.isQ[j,n]:
				p = 0
				for k in range(case.nBus):
					if case.isT[k]:
						if j==n:
							if j==k: M[i,p] = V[j]*sum([LambdaP(V[m],case.G[j,m],case.B[j,m],theta[j] - theta[m]) for m in case.K[j]])
							elif (k in case.K[j] and k!=j): M[i,p] = -V[j]*LambdaP(V[k],case.G[j,k],case.B[j,k],theta[j] - theta[k])
						else:
							if k == j: M[i,p] = V[j]*LambdaP(V[n],self.G[j,n],self.B[j,n],theta[j] - theta[n])
							else: M[i,p] = -V[j]*LambdaP(V[n],self.G[j,n],self.B[j,n],theta[j] - theta[n])
						#print(' --> dQdT[{},{}] = dQ({},{})/dT({}) = {}'.format(i,p,j,n,k,M[i,p]))
						p+=1				
				i+=1
	#print(M)
	return M

def dQdV(V,theta,case): #{{{2
	L = np.zeros((int(case.isQ.sum()), int(case.isV.sum())))
	i = 0
	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isQ[j,n]:
				p = 0
				for k in range(case.nBus):
					if case.isV[k]:
						if j==n:
							if j==k: L[i,p] = -2*V[j]*case.B[j,j] + sum([LambdaP(V[m],case.G[j,m],case.B[j,m],theta[j] - theta[m]) for m in case.K[j]])
							elif (k in case.K[j]): L[i,p] = LambdaP(V[j],-case.B[j,k],case.G[j,k],theta[j] - theta[k])
						else: # Calculating dQjn/dVk
							if k == j: L[i,p] = -2*V[j]*case.Bsh[j,n] + LambdaP(V[n],-self.B[j,n],self.G[j,n],theta[j] - theta[n])
							else: L[i,p] = LambdaP(V[n],-self.B[j,n],self.G[j,n],theta[j] - theta[n])

						#print(' --> dQdV[{},{}] = dQ({},{})/dV({}) = {}'.format(i,p,j,n,k,L[i,p]))
						p+=1
				i+=1
	#print(L)
	return L

def h(V,theta,case): #{{{2
	h = np.zeros((int(case.isP.sum()) + int(case.isQ.sum()),1))
	i = 0
	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isP[j,n]:
				#print(' >>>>>>>>>> h({}) = P({},{})'.format(i,j,n))
				if j==n: # If measuring power injection on bar j
					h[i] = V[j]**2*case.G[j,j] + V[j]*sum([LambdaP(V[m],case.G[j,m],case.B[j,m],theta[j] - theta[m]) for m in case.K[j]])
				else:	# If measuring power flux from bar j to n
					h[i] = V[j]**2*case.Gsh[j,n] + V[j]*LambdaP(V[n],case.G[j,n],case.B[j,n],theta[j] - theta[n])
				i +=1
	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isQ[j,n]:
				#print(' >>>>>>>>>> h({}) = Q({},{})'.format(i,j,n))
				if j==n:
					h[i] = -V[j]**2*case.B[j,j] + V[j]*sum([ LambdaP(V[m],-case.B[j,m],case.G[j,m],theta[j] - theta[m]) for m in case.K[j]])
				else:
					h[i] = -V[j]**2*case.Bsh[j,n] + V[j]*LambdaP(V[m],-case.B[j,m],case.G[j,m],theta[j] - theta[m])
				i +=1
	return h

def Z(P,Q,case): #{{{2
	i=0
	Z = np.empty((int(case.isP.sum()) + int(case.isQ.sum()),1))
	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isP[j,n]:
				Z[i] = P[j,n]
				#print(' >>>>>>>>>> Z({}) = P({},{}) = {}'.format(i,j,n,Z[i]))
				i+=1

	for j in range(case.nBus):	
		for n in range(case.nBus):
			if case.isQ[j,n]:
				Z[i] = Q[j,n]
				#print(' >>>>>>>>>> Z({}) = Q({},{}) = {}'.format(i,j,n,Z[i]))
				i += 1
	return Z

def Jac(V,theta,case): #{{{2

	H = dPdT(V,theta,case)
	N = dPdV(V,theta,case)
	M = dQdT(V,theta,case)
	L = dQdV(V,theta,case)

	dP = conc((H,N),axis=1)
	dQ = conc((M,L),axis=1)

	return conc((dP,dQ),axis=0)

#}}}1

# Functions for Droop-specific power flow {{{1 
