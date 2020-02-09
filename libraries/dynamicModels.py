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
def SM2(k,x,I,genData):
	#F = np.zeros(2)
	omega = x[k]
	delta = x[k+1]
	H = genData[k,2]
	D = genData[k,3]
	pm0 = genData[k,17]
	Elq0 = genData[k,19]
	Iq = np.real(I)
	F[k] = ( pm0 - Elq0*Iq - D*omega )/(2*H)	# omega = x[k+1]
	F[k+1] = omega
	return F

# Synchronous machine one-axis (third-order) model
def SM1A:
	EPq = x[k]
	omega = x[k+1]
	delta = x[k+2]
	EFD0 = genData[k,18]
	tPdo = genData[k,10]
	F[k] = 1/tPdo*(EFD0 - EPq + (xd - xPd)*Id)
	F[k+1] = ( pm0 - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	F[k+2] = x[k+1]
