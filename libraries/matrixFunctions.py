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

def dPdT(V,theta,K,a,y,Y,bsh,isP): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	nbus = len(isP)
	H = np.empty((numP, nbus))	# Preallocating
	i = 0	# Starting counter
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				for k in range(nbus):
					if j==n: # If measuring power injection
						if k==j: # If calculating on angle j
							H[i,k] = V[j]*sum([V[m]*(-G[j,m]*sin(theta[j] - theta[m]) + B[j,m]*cos(theta[j] - theta[m])) for m in K[j]])
						elif (k in K[j]): # If calculating on angle k, n and j are connected
							H[i,k] = V[j]*V[k]*(G[j,k]*sin(theta[j] - theta[k]) - B[j,k]*cos(theta[j] - theta[k]))
						else: # If n and j are not connected
							H[i,k] = 0
					else:	# If measuring power flow
						if k!=j and k!=n: # Since power transfer is function of only angles j and n
							H[i,k] = 0
						else:
							if k==j: H[i,k] = V[j]*V[n]*a[j,n]*a[n,j]*(  g[j,n]*sin(theta[j] - theta[n]) - b[j,n]*cos(theta[j] - theta[n]))
							else: H[i,k] =    V[j]*V[n]*a[j,n]*a[n,j]*( -g[j,n]*sin(theta[j] - theta[n]) + b[j,n]*cos(theta[j] - theta[n]))
							
				i+=1
	# Deleting the reference bar
	H = np.delete(H,0,axis=1)
	return H

def dPdV(V,theta,K,a,y,Y,bsh,isP):#{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	nbus = len(isP)
	N = np.empty((numP, nbus))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				for k in range(nbus):
					if j==n:
						if k==j: 
							N[i,k] = 2*V[j]*G[j,j] + sum([V[m]*(G[j,m]*cos(theta[j] - theta[m]) + B[j,m]*sin(theta[j] - theta[m])) for m in K[j]])
						elif k in K[j]: 
							N[i,k] = V[j]*(G[j,k]*cos(theta[j] - theta[k]) + B[j,k]*sin(theta[j] - theta[k]))
						else: N[i,k] = 0
					else:
						if k!=j and k!=n:
							N[i,k] = 0
						else:
							if k == j: N[i,k] = 2*a[j,n]*V[j]*g[j,n] - a[j,n]*a[n,j]*V[n]*(g[j,n]*cos(theta[j] - theta[n]) + b[j,n]*sin(theta[j] - theta[n]))
							else: N[i,k] = -a[j,n]*a[n,j]*V[j]*(g[j,n]*cos(theta[j] - theta[n]) + b[j,n]*sin(theta[j] - theta[n]))
				i+=1
	return N

def dQdT(V,theta,K,a,y,Y,bsh,isP): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	nbus = len(isP)
	M = np.empty((numP, nbus))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				for k in range(nbus):
					if j==n: # If measuring power injection
						if j==k: # If calculating on angle j
							M[i,k] = V[j]*sum([V[m]*(G[j,m]*cos(theta[j] - theta[m]) + B[j,m]*sin(theta[j] - theta[m])) for m in K[j]])
						elif (k in K[j] and k!=j):
							M[i,k] = -V[j]*V[k]*(G[j,k]*cos(theta[j] - theta[k]) + B[j,k]*sin(theta[j] - theta[k]))
						else: M[i,k] = 0
					else:
						if k!=j and k!=n: #Since power transfer is function of only angles j and n
							M[i,k] = 0
						else:
							if k == j: M[i,k] = -a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n]) + b[j,n]*sin(theta[j] - theta[n]))
							else: M[i,k] = a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n]) + b[j,n]*sin(theta[j] - theta[n]))
							
				i+=1
	M = np.delete(M,0,axis=1)
	return M

def dQdV(V,theta,K,a,y,Y,bsh,isP): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	nbus = len(isP)
	L = np.ones((numP, nbus))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				for k in range(nbus):
					if j==n: # If measuring power injection
						if j==k: L[i,k] = -2*V[j]*B[j,j] + sum([V[m]*(G[j,m]*sin(theta[j] - theta[m]) - B[j,m]*cos(theta[j] - theta[m])) for m in K[j]])
						elif (k in K[j]): L[i,k] = V[j]*(G[j,k]*sin(theta[j] - theta[k]) - B[j,k]*cos(theta[j] - theta[k]))
						else: L[i,k] = 0
					else:
						if k!=j and k!=n: #Since power transfer is function of only voltages j and n
							L[i,k] = 0
						else:
							if k == j: L[i,k] = -2*a[j,n]**2*V[j]*(b[j,n] + bsh[j,n]) + a[j,n]*a[n,j]*V[n]*(b[j,n]*cos(theta[j] - theta[n]) - g[j,n]*sin(theta[j] - theta[n]))
							else: L[i,k] = a[j,n]*a[n,j]*V[j]*(b[j,n]*cos(theta[j] - theta[n]) - g[j,n]*sin(theta[j] - theta[n]))
				i += 1
	return L

def h(V,theta,K,a,y,Y,bsh,isP,isV): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	numV = int(isV.sum())
	nbus = len(isP)
	h = np.zeros((2*numP + numV,1))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				if j==n: # If measuring power injection on bar j
					h[i] = V[j]**2*G[j,j] + V[j]*sum([ V[m]*(G[j,m]*cos(theta[j] - theta[m]) + B[j,m]*sin(theta[j] - theta[m])) for m in K[j]])
				else:	# If measuring power flux from bar j to n
					h[i] = V[j]**2*a[j,n]**2*g[j,n] - a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n]) + b[j,n]*sin(theta[j] - theta[n]))
					
				i +=1
				
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				if j==n:
					h[i] = -V[j]**2*bsh[j,j] + sum([ -V[j]**2*a[j,m]**2*(b[j,m] + bsh[j,m]) + a[m,j]*a[j,m]*V[j]*V[m]*( -g[j,m]*sin(theta[j] - theta[m]) + b[j,m]*cos(theta[j] - theta[m]) ) for m in K[j]])
				else:
					h[i] =  -V[j]**2*a[j,n]**2*(b[j,n] + bsh[j,n]) + a[n,j]*a[j,n]*V[j]*V[n]*( -g[j,n]*sin(theta[j] - theta[n]) + b[j,n]*cos(theta[j] - theta[n]) )
				i +=1

	for j in range(nbus):	
		if isV[j]:
			h[i] = V[j]
			i += 1
	return h

def Z(P,Q,isP,V,isV): #{{{1
	numP = int(isP.sum())
	numV = int(isV.sum())
	nbus = int(len(V))
	i=0
	Z = np.empty((2*numP + numV,1))
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				Z[i] = P[j,n]
				i+=1
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				Z[i] = Q[j,n]
				i += 1

	for j in range(nbus):	
		if isV[j]:
			Z[i] = V[j]
			i += 1

	return Z

def Jac(V,theta,K,a,y,Y,bsh,isP,isV): #{{{1

	nbus = len(isP)
	numP = isP.sum()
	numV = isV.sum()

	H = dPdT(V,theta,K,a,y,Y,bsh,isP)
	N = dPdV(V,theta,K,a,y,Y,bsh,isP)
	M = dQdT(V,theta,K,a,y,Y,bsh,isP)
	L = dQdV(V,theta,K,a,y,Y,bsh,isP)

	O = np.identity(nbus)
	deleteList = array([i for i in range(nbus) if isV[i]==0])
	O = np.delete(O,deleteList,axis=0) # Deleting non-V measures
	O = conc((np.zeros((O.shape[0],O.shape[1]-1)),O),axis=1)

	dP = conc(( H, N),axis=1)
	dQ = conc(( M, L),axis=1)
	dPdQ = conc((dP,dQ),axis=0)
	
	return conc((dPdQ,O),axis=0)
