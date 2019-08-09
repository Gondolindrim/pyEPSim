# -------------------------------------------------
# UNIVERSITY OF SAO PAULO
# SÃO CARLOS SCHOOL OF ENGINEERING (EESC)
# DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING (SEL)
# TITLE: pyEPSim dynamic models library
# AUTHOR: Álvaro Augusto "Gondolindrim" Volpato
# DATE: 04/07/2018
# VERSION: 1.3
# DESCRIPTION: this file contains the dynamic models of the synchronous generators used for the dynamic simulation method of pyEPSim.
# -------------------------------------------------

# Synchronous machine second-order model with damping
def sm2d(x,t,C,D,Yred,V,pm,genData):
	nGen = genData.shape[0]
	F = np.zeros(2*nGen)
	for k in range(nGen):
		F[2*k] = x[2*k+1]	# delta} = x[k]
		F[2*k+1] = ( pm[k] - (V[k]**2)*real(Yred[k,k]) - sum( [C[k,j]*sin(x[2*k] - x[2*j]) + D[k,j]*cos(x[2*k] - x[2*j]) for j in range(nGen)]) - genData[k,3]*x[2*k+1] )/(2*genData[k,2])	# omega = x[k+1]

	return F

# Synchronous machine second-order model without damping
def sm2(x,t,C,D,Yred,V,pm,genData):
	nGen = genData.shape[0]
	F = np.zeros(2*nGen)
	for k in range(nGen):
		F[2*k] = x[2*k+1]	# delta} = x[k]
		F[2*k+1] = ( pm[k] - (V[k]**2)*real(Yred[k,k]) - sum( [C[k,j]*sin(x[2*k] - x[2*j]) + D[k,j]*cos(x[2*k] - x[2*j]) for j in range(nGen)]) )/(2*genData[k,2])	# omega = x[k+1]

	return F
