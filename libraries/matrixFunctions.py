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

def dPdT(V,theta,K,a,phi,y,Y,bsh,isP,isT): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	numT = int(isT.sum())
	nbus = len(isP)
	H = np.empty((numP, numT))	# Preallocating
	i = 0	# Starting counter
	for j in range(nbus):
		for n in range(nbus):
			if isP[j,n]:
				p = 0
				for k in range(nbus):
					if isT[k]:
						if j==n: # If measuring power injection
							if k==j: # If calculating on angle j
								H[i,p] = V[j]*sum([V[m]*(-G[j,m]*sin(theta[j] - theta[m]) + B[j,m]*cos(theta[j] - theta[m])) for m in K[j]])
							elif (k in K[j]): # If calculating on angle k, n and j are connected
								H[i,p] = V[j]*V[k]*(G[j,k]*sin(theta[j] - theta[k]) - B[j,k]*cos(theta[j] - theta[k]))
							else: # If n and j are not connected
								H[i,p] = 0
						else:	# If measuring power flow
							if k!=j and k!=n: # Since power transfer is function of only angles j and n
								H[i,p] = 0
							else:
								if k==j: H[i,p] = V[j]*V[n]*a[j,n]*a[n,j]*(  g[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]) - b[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
								else: H[i,p] =    V[j]*V[n]*a[j,n]*a[n,j]*( -g[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
						#print(' --> dPdT[{},{}] = dP({},{})/dT({})'.format(i,k,j,n,k))
						p+=1
				i+=1
	# Deleting the reference bar
	#print(H)
	return H

def dPdV(V,theta,K,a,phi,y,Y,bsh,isP,isV):#{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	numV = int(isV.sum())
	nbus = len(isP)
	N = np.empty((numP, numV))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				p = 0
				for k in range(nbus):
					if isV[k]:
						if j==n:
							if k==j: N[i,p] = 2*V[j]*G[j,j] + sum([V[m]*(G[j,m]*cos(theta[j] - theta[m]) + B[j,m]*sin(theta[j] - theta[m])) for m in K[j]])
							elif k in K[j]: N[i,p] = V[j]*(G[j,k]*cos(theta[j] - theta[k]) + B[j,k]*sin(theta[j] - theta[k]))
							else: N[i,p] = 0
						else:
							if k!=j and k!=n: N[i,p] = 0
							else:
								if k == j: N[i,p] = 2*a[j,n]**2*V[j]*g[j,n] - a[j,n]*a[n,j]*V[n]*(g[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
								else: N[i,p] = -a[j,n]*a[n,j]*V[j]*(g[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
						#print(' --> dPdV[{},{}] = dP({},{})/dV({}) = {}'.format(i,p,j,n,k,N[i,p]))
						p+=1
				i+=1
	#print(N)
	return N

def dQdT(V,theta,K,a,phi,y,Y,bsh,isQ,isT): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numQ = int(isQ.sum())
	numT = int(isT.sum())
	nbus = len(isQ)
	M = np.empty((numQ, numT))
	i = 0
	for j in range(nbus):
		for n in range(nbus):
			if isQ[j,n]:
				p = 0
				for k in range(nbus):
					if isT[k]:
						if j==n: # If measuring power injection
							if j==k: M[i,p] = V[j]*sum([V[m]*(G[j,m]*cos(theta[j] - theta[m] + phi[j,n] - phi[n,j]) + B[j,m]*sin(theta[j] - theta[m] + phi[j,n] - phi[n,j])) for m in K[j]])
							elif (k in K[j] and k!=j): M[i,p] = V[j]*V[k]*(-G[j,k]*cos(theta[j] - theta[k]) - B[j,k]*sin(theta[j] - theta[k]))
							else: M[i,p] = 0
						else:
							if k!=j and k!=n: M[i,p] = 0
							else:
								if k == j: M[i,p] = -a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
								else: M[i,p] = a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
						#print(' --> dQdT[{},{}] = dQ({},{})/dT({}) = {}'.format(i,p,j,n,k,M[i,p]))
						p+=1				
				i+=1
	#print(M)
	return M

def dQdV(V,theta,K,a,phi,y,Y,bsh,isQ,isV): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isQ.sum())
	numV = int(isV.sum())
	nbus = len(isQ)
	L = np.ones((numP, numV))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isQ[j,n]:
				p = 0
				for k in range(nbus):
					if isV[k]:
						if j==n:
							if j==k: L[i,p] = -2*V[j]*B[j,j] + sum([V[m]*(G[j,m]*sin(theta[j] - theta[m]) - B[j,m]*cos(theta[j] - theta[m])) for m in K[j]])
							elif (k in K[j]): L[i,p] = V[j]*(G[j,k]*sin(theta[j] - theta[k]) - B[j,k]*cos(theta[j] - theta[k]))
							else: L[i,p] = 0
						else:
							if k!=j and k!=n: L[i,p] = 0
							else:
								if k == j: L[i,p] = -2*a[j,n]**2*V[j]*(b[j,n] + bsh[j,n]) + a[j,n]*a[n,j]*V[n]*(b[j,n]*cos(theta[j] - theta[n]) - g[j,n]*sin(theta[j] - theta[n]))
								else: L[i,p] = a[j,n]*a[n,j]*V[j]*(b[j,n]*cos(theta[j] - theta[n]) - g[j,n]*sin(theta[j] - theta[n]))
						#print(' --> dQdV[{},{}] = dQ({},{})/dV({}) = {}'.format(i,p,j,n,k,L[i,p]))
						p+=1
				i+=1
	#print(L)
	return L

def h(V,theta,K,a,phi,y,Y,bsh,isP,isQ,isV): #{{{1
	G = real(copy(Y))
	B = imag(copy(Y))
	g = real(copy(y))
	b = imag(copy(y))
	numP = int(isP.sum())
	numQ = int(isQ.sum())
	nbus = len(isP)
	h = np.zeros((numP + numQ,1))
	i = 0
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				#print(' >>>>>>>>>> h({}) = P({},{})'.format(i,j,n))
				if j==n: # If measuring power injection on bar j
					h[i] = V[j]**2*G[j,j] + V[j]*sum([ V[m]*(G[j,m]*cos(theta[j] - theta[m]) + B[j,m]*sin(theta[j] - theta[m])) for m in K[j]])
				else:	# If measuring power flux from bar j to n
					h[i] = (a[j,n]*V[j])**2*g[j,n] - a[j,n]*a[n,j]*V[j]*V[n]*(g[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]))
					
				i +=1
	for j in range(nbus):	
		for n in range(nbus):
			if isQ[j,n]:
				#print(' >>>>>>>>>> h({}) = Q({},{})'.format(i,j,n))
				if j==n:
					h[i] = -V[j]**2*B[j,j] + V[j]*sum([ V[m]*(G[j,m]*sin(theta[j] - theta[m]) - B[j,m]*cos(theta[j] - theta[m])) for m in K[j]])
				else:
					h[i] = -(a[j,n]*V[j])**2**(b[j,n] + bsh[j,n]) + a[j,n]*a[n,j]*V[j]*V[n]*( -g[j,n]*sin(theta[j] - theta[n] + phi[j,n] - phi[n,j]) + b[j,n]*cos(theta[j] - theta[n] + phi[j,n] - phi[n,j]) )
				i +=1
	return h

def Z(P,Q,isP,isQ): #{{{1
	numP = int(isP.sum())
	numQ = int(isQ.sum())
	nbus = int(len(isP))
	i=0
	Z = np.empty((numP + numQ,1))
	for j in range(nbus):	
		for n in range(nbus):
			if isP[j,n]:
				Z[i] = P[j,n]
				#print(' >>>>>>>>>> Z({}) = P({},{}) = {}'.format(i,j,n,Z[i]))
				i+=1

	for j in range(nbus):	
		for n in range(nbus):
			if isQ[j,n]:
				Z[i] = Q[j,n]
				#print(' >>>>>>>>>> Z({}) = Q({},{}) = {}'.format(i,j,n,Z[i]))
				i += 1
	return Z

def Jac(V,theta,K,a,phi,y,Y,bsh,isP,isQ,isV,isT): #{{{1

	numP = isP.sum()
	numQ = isQ.sum()
	numV = isV.sum()
	numT = isT.sum()

	nbus = len(isP)

	H = dPdT(V,theta,K,a,phi,y,Y,bsh,isP,isT)
	N = dPdV(V,theta,K,a,phi,y,Y,bsh,isP,isV)
	M = dQdT(V,theta,K,a,phi,y,Y,bsh,isQ,isT)
	L = dQdV(V,theta,K,a,phi,y,Y,bsh,isQ,isV)

	dP = conc(( H, N),axis=1)
	dQ = conc(( M, L),axis=1)
	dPdQ = conc((dP,dQ),axis=0)
	
	return dPdQ
