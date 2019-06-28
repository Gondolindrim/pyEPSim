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

def dPdT(V,theta,K,a,y,Y,bsh,isP): #{{{2
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

def dPdV(V,theta,K,a,y,Y,bsh,isP):#{{{2
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

def dQdT(V,theta,K,a,y,Y,bsh,isP): #{{{2
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

def dQdV(V,theta,K,a,y,Y,bsh,isP): #{{{2
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

def h(V,theta,K,a,y,Y,bsh,isP,isV): #{{{2
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

def Z(P,Q,isP,V,isV): #{{{2
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

def Jac(V,theta,K,a,y,Y,bsh,isP,isV): #{{{2

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

def reduceGrid(Y,Yload,V,genData): #{{{2
	nGen = genData.shape[0]
	nBus = Y.shape[0]

	# Building component matrixes

	Y1 = Y[0:nGen,0:nGen]
	Y2 = Y[0:nGen,nGen:nBus+1]
	Y3 = Y[nGen:nBus,0:nGen]
	Y4 = Y[nGen:nBus+1,nGen:nBus+1]

	Ylg = np.diag(Yload[0:nGen])

	Yll = np.diag(Yload[nGen:nBus])

	Ytrans = np.array([1/(data[4] + 1j*data[7]) for data in genData]) # Y' no livro
	Ytrans = np.diag(Ytrans)
	
	YA = Ytrans

	YB = conc((-Ytrans,np.zeros((nGen,nBus-nGen))),axis=1)
	
	YC = conc((-Ytrans,np.zeros((nBus-nGen,nGen))),axis=0)

	YDtop = conc((Ytrans+Y1+Ylg,Y2),axis=1)
	YDbot = conc((Y3,Y4 + Yll),axis=1)

	YD = conc((YDtop,YDbot),axis=0)

	# Calculating YRED and C and D coefficients

	Yred = YA - YB @ inv(YD) @ YC

	C = np.zeros((nGen,nGen))
	D = np.zeros((nGen,nGen))
	
	for i in range(nGen):
		for j in range(nGen):
			if j != i:
				C[i,j] += V[i]*V[j]*imag(Yred[i,j])
				D[i,j] += V[i]*V[j]*real(Yred[i,j])

	return [Yred,C,D]

def odeFault(x,t,C,D,Yred,V,pm,genData):
	nGen = genData.shape[0]
	F = np.zeros(2*nGen)
	for k in range(nGen):
		F[2*k] = x[2*k+1]	# delta} = x[k]
		F[2*k+1] = ( pm[k] - (V[k]**2)*real(Yred[k,k]) - sum( [C[k,j]*sin(x[2*k] - x[2*j]) + D[k,j]*cos(x[2*k] - x[2*j]) for j in range(nGen)])  )/(2*genData[k,2]) - genData[k,3]*x[2*k+1]/(2*genData[k,2])	# omega = x[k+1]

	return F
