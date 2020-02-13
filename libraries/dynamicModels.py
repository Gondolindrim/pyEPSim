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
import numpy as np
# Algebraic function for solving V and I
def solveI(x,*args): #{{{
	Elq, Eld, ra, xPq, xPd, Y, k, V = args
	Iq = x[0]
	Id = x[1]
	Vq = x[2]
	Vd = x[3]

	V[k] = Vq + 1j*Vd

	F0 = Iq - np.real(sum([Y[k,m]*V[m] for m in range(len(Y))]))
	F1 = Id - np.imag(sum([Y[k,m]*V[m] for m in range(len(Y))]))
	F2 = -Vq + Elq - ra*Iq + xPd*Id
	F3 = -Vd + Eld - ra*Id - xPq*Iq
	return [F0, F1, F2, F3] #}}}

# Synchronous machine second-order model with damping
def SMC(x,*args): #{{{1
	Ig, genData = args

	omega = x[0]
	delta = x[1]
	
	Iq = np.real(Ig)
	Id = np.imag(Ig)

	xPd = genData.xPd
	xPq = genData.xPq
	H = genData.H
	D = genData.D
	r = genData.ra
	pm = genData.pm0	# Mechanical power is constant equal to its initial value
	el0 = genData.el0
	Elq = np.real(genData.el0)	# Same as mechanical power
	Eld = np.imag(genData.el0)
	
	omegaP = ( pm - Elq*Iq - D*omega )/(2*H)	# omega = x[k+1]	
	deltaP = omega
	return [omegaP, deltaP] #}}}1

# Synchronous machine one-axis (third-order) model
def SM1A(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)

	Elq = x[0]
	Eld = 0
	omega = x[1]
	delta = x[2]

	H = genData.H
	D = genData.D
	EFD0 = genData.efd0
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	pm = genData.pm0	# Mechanical power is constant equal to its initial value

	dElq = 1/tPdo*(EFD0 - Elq - (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta = omega

	return [dElq, domega, ddelta] #}}}1

# Synchronous machine one-axis (third-order) model with turbine and governor equations
def SM1A_TUR_GOV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	Elq = x[0]
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	
	H = genData.H
	D = genData.D
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	EFD0 = genData.efd0
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	tG = genData.tG
	tT = genData.tT

	dElq = 1/tPdo*(EFD0 - Elq - (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)

	return [dElq, domega, ddelta, dpm, dpPm] #}}}1

# Synchronous machine one-axis (third-order) model with turbine, governor, AVR and PSS equations
def SM1A_TUR_GOV_AVR_PSS(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = np.imag(genData.el0)
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	vAVR = x[5]
	vWash = x[6]
	vPSS = x[7]
	
	H = genData.H
	D = genData.D
	r = genData.ra
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	EFD0 = genData.efd0
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	tG = genData.tG
	tT = genData.tT
	Vt0 = genData.vRef

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	#print('\n >>> Gen {} passed EL = {}'.format(genData.busNumber, Elq + 1j*Eld))
	#print(' >>> Gen {} passed I = {}'.format(genData.busNumber, Iq + 1j*Id))
	#print(' >>> Gen {} model Vt = {}'.format(genData.busNumber,Vq + 1j*Vd))
	#print(' >>> Gen {} model Vt0 = {}'.format(genData.busNumber,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	EFD = EFD0 + vAVR + vPSS

	dElq = 1/tPdo*(EFD - Elq - (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dvAVR = (Ke*(Vt - Vt0) - vAVR)/Te

	#print(' >>> Gen {} model Vt - Vtref = {}\n'.format(genData.busNumber,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	k = genData.busNumber

	return [dElq, domega, ddelta, dpm, dpPm, dvAVR, dvWash, dvPSS] #}}}1

# Synchronous machine two-axis (fourth-order) model with turbine, governor, AVR and PSS equations
def SM2A_TUR_GOV_AVR_PSS(x,*args): #{{{1
	k, V, Y, genData = args
	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]
	vAVR = x[6]
	vWash = x[7]
	vPSS = x[8]
	
	H = genData.H
	D = genData.D
	r = genData.ra
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	EFD0 = genData.efd0
	tPdo = genDatatPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT
	Vt0 = genData.vRef
	Vt = V[k]
	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	EFD = EFD0 + vAVR + vPSS

	Iq, Id, Vq, Vd = fsolve(solveI, [1, 0, 1, 0], args = (Elq, Eld, xPq, xPd, Y, k, V))

	dEld = 1/tPqo*(Eld + (xq - xPq)*Iq)
	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpmSet = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dvAVR = (Ke*(Vt - Vt0) - vAVR)/Te
	dvWash = KPss/(2*H)*( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega ) - vWash/Tw
	dvPSS = T2/T1*( KPss/(2*H)*( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega ) +  vWash*(1/T2 - 1/Tw) - 1/T2*vPSS)

	return [dElq, dEld, domega, ddelta, dpm, dpmSet, dvAVR, dvWash, dvPSS] #}}}1
