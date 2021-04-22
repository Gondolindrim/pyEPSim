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
def solveCurrents(y, *args):	#{{{1
	Y, EL, angles, case = args
	#print('\n >>> Passed angleReferences = {}'.format(angles))
	#print(' >>> Passed EL = {}'.format(EL))

	nGen  = case.nGen

	F = np.zeros(2*nGen, dtype = float)
	V = [y[k] + 1j*y[k + nGen] for k in range(nGen)]

	I = Y @ V
	#print(I)
	for k in range(case.nGen):
		# I and V are the current and voltage vectors in the synchronous reference
		gen = case.genData[k]
		ELq = np.real(EL[k])
		ELd = np.imag(EL[k])
	
		Imachine = I[k]*np.e**(-1j*angles[k])
		Vmachine = V[k]*np.e**(-1j*angles[k])

		#print(' >>> Imachine = {}'.format(Imachine))
		Iq = np.real(Imachine)
		Id = np.imag(Imachine)

		Vq = np.real(Vmachine)
		Vd = np.imag(Vmachine)
		#print(Vq + 1j*Vd)
		
		F[k] = Vq - ELq + gen.ra*Iq - gen.xPd*Id
		F[k + nGen] = Vd - ELd + gen.ra*Id + gen.xPq*Iq
		
	#print('>>> Calculated F function norm = {}\n'.format(F))
	#input()
	return F	#}}}

# MODEL (1): synchronous machine second-order model with damping
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
	pm = genData.pm0	# Mechanical power is constant equal to its initial value
	el0 = genData.el0
	Elq = np.real(genData.el0)	# Same as mechanical power
	Eld = np.imag(genData.el0)
	
	domega = ( pm - Elq*Iq - ELd*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)	
	ddelta = omega
	return [domega, ddelta] #}}}1

# MODEL (2): synchronous machine one-axis (third-order) model
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

	dElq = 1/tPdo*(EFD0 - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - ELd*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta = omega

	return [dElq, domega, ddelta] #}}}1

# MODEL (3): synchronous machine one-axis (third-order) model with turbine and governor equations
def SM1A_TUR_GOV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	Elq = x[0]
	Eld = np.imag(genData.el0)
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

	dElq = 1/tPdo*(EFD0 - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)

	return [dElq, domega, ddelta, dpm, dpPm] #}}}1

# MODEL (4): synchronous machine one-axis (third-order) model with turbine, governor, AVR and PSS equations
def SM1A_TUR_GOV_AVR_PSS(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = np.imag(genData.el0)
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	EFD = x[5]
	vWash = x[6]
	vPSS = x[7]
	
	H = genData.H
	D = genData.D
	r = genData.ra
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	tG = genData.tG
	tT = genData.tT
	Vtref = genData.vRef

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vt - Vtref + vPSS) - EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	return [dElq, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS] #}}}1

# MODEL (5): synchronous machine one-axis (third-order) model with turbine, governor, AVR and PSS equations, and PF-QV droop control
def SM1A_TUR_GOV_AVR_PSS_PFQV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = np.imag(genData.el0)
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	EFD = x[5]
	vWash = x[6]
	vPSS = x[7]
	
	H = genData.H
	D = genData.D
	r = genData.ra
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	tG = genData.tG
	tT = genData.tT

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	# Droop equations
	P = Vq*Iq + Vd*Id
	Q = Vd*Iq - Vq*Id
	pmSet = genData.P0 - genData.kP * omega
	Vtref = genData.vRef - (Q - genData.Q0)/genData.kQ

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) - EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	return [dElq, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS] #}}}1

# MODEL (6): synchronous machine one-axis (third-order) model with turbine, governor, AVR and PSS equations, and PF-Q(dV) droop control
def SM1A_TUR_GOV_AVR_PSS_PFQdV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = np.imag(genData.el0)
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	EFD = x[5]
	vWash = x[6]
	vPSS = x[7]
	Vtref = x[8]	

	H = genData.H
	D = genData.D
	r = genData.ra
	tPdo = genData.tPdo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	tG = genData.tG
	tT = genData.tT

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	# Droop equations
	P = Vq*Iq + Vd*Id
	Q = Vd*Iq - Vq*Id
	pmSet = genData.P0 - genData.kP * omega

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) - EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	dVtref = - (Q - genData.Q0)/genData.kQ

	return [dElq, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS, dVtref] #}}}1

# MODEL (7): synchronous machine two-axis (fourth-order) model with turbine, governor, AVR and PSS equations
def SM2A_TUR_GOV_AVR_PSS(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)

	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]
	EFD = x[6]
	vWash = x[7]
	vPSS = x[8]
	
	H = genData.H
	D = genData.D
	r = genData.ra
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	tPdo = genData.tPdo
	tPqo = genData.tPqo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT
	Vt0 = genData.vRef
	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2
	Ke = genData.Ke
	Te = genData.Te
	KPss = genData.KPss

	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	dEld = -1/tPqo*(Eld + (xq - xPq)*Iq)
	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) - EFD)/Te
	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	return [dElq, dEld, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS] #}}}1

#MODEL (8): synchronous machine two-axis (fourth-order) model with turbine, governor, AVR and PSS equations, and PF-QV droop control
def SM2A_TUR_GOV_AVR_PSS_PFQV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]
	EFD = x[6]
	vWash = x[7]
	vPSS = x[8]

	H = genData.H
	D = genData.D
	r = genData.ra
	tPdo = genData.tPdo
	tPqo = genData.tPqo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	# Droop equations
	P = Vq*Iq + Vd*Id
	Q = Vd*Iq - Vq*Id
	pmSet = genData.P0 - genData.kP * omega
	Vtref = genData.vRef - (Q - genData.Q0)/genData.kQ

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	dEld = -1/tPqo*(Eld + (xq - xPq)*Iq)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) - EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	return [dElq, dEld, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS] #}}}1

#MODEL (9): synchronous machine two-axis (fourth-order) model with turbine, governor, AVR and PSS equations, and PF-Q(dV) droop control
def SM2A_TUR_GOV_AVR_PSS_PFQdV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]
	EFD = x[6]
	vWash = x[7]
	vPSS = x[8]
	Vtref = x[9]
	Q0 = x[10]

	H = genData.H
	D = genData.D
	r = genData.ra
	tPdo = genData.tPdo
	tPqo = genData.tPqo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	# Droop equations
	P = Vq*Iq + Vd*Id
	Q = Vd*Iq - Vq*Id
	pmSet = genData.P0 - genData.kP * omega

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	dEld = -1/tPqo*(Eld + (xq - xPq)*Iq)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) + EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	dVtref =  -(Q - Q0)/genData.kQ
	dQ0 = -genData.kReg*genData.ratedPower*(dVt0)
	return [dElq, dEld, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS, dVt0, dQ0] #}}}1

# MODEL (10): synchronous machine two-axis (fourth-order) model with turbine and governor equations
def SM2A_TUR_GOV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)

	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]	
	
	H = genData.H
	D = genData.D
	r = genData.ra
	pmSet = genData.pm0	# Initial value of Pm setpoint is the initial power
	tPdo = genData.tPdo
	tPqo = genData.tPqo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT

	EFD = genData.efd0

	dEld = -1/tPqo*(Eld + (xq - xPq)*Iq)
	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - Id*Eld - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	
	return [dElq, dEld, domega, ddelta, dpm, dpPm] #}}}1

#MODEL (11): synchronous machine two-axis (fourth-order) model with turbine, governor, AVR and PSS equations, and iPF-Q(dV) droop control
def SM2A_TUR_GOV_AVR_PSS_iPFQdV(x,*args): #{{{1
	Ig, genData = args
	Iq, Id = np.real(Ig), np.imag(Ig)
	
	Elq = x[0]
	Eld = x[1]
	omega = x[2]
	delta = x[3]
	pm = x[4]
	pPm = x[5]
	EFD = x[6]
	vWash = x[7]
	vPSS = x[8]
	Vtref = x[9]
	Q0 = x[10]
	P0 = x[11]

	H = genData.H
	D = genData.D
	r = genData.ra
	tPdo = genData.tPdo
	tPqo = genData.tPqo
	xPd = genData.xPd
	xPq = genData.xPq
	xd = genData.xd
	xq = genData.xq
	tG = genData.tG
	tT = genData.tT

	KPss = genData.KPss
	Ke = genData.Ke
	Te = genData.Te
	
	Vq = Elq - r*Iq + xPd*Id
	Vd = Eld - r*Id - xPq*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)

	# Droop equations
	P = Vq*Iq + Vd*Id
	Q = Vd*Iq - Vq*Id
	pmSet = P0 - genData.kP * omega

	#print('\n >>> Gen \'{}\' passed EL = {}'.format(genData.busName, Elq + 1j*Eld))
	#print(' >>> Gen \'{}\' passed I = {}'.format(genData.busName, Iq + 1j*Id))
	#print(' >>> Gen \'{}\' model Vt = {}'.format(genData.busName,Vq + 1j*Vd))
	#print(' >>> Gen \'{}\' model Vt0 = {}'.format(genData.busName,Vt0))

	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	dElq = 1/tPdo*(EFD - Elq + (xd - xPd)*Id)
	dEld = -1/tPqo*(Eld + (xq - xPq)*Iq)
	domega = ( pm - Elq*Iq - Eld*Id - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dEFD = (Ke*(Vtref - Vt + vPSS) - EFD)/Te

	#print(' >>> Gen at bus \'{}\' model Vt - Vtref = {}\n'.format(genData.busName,Vt - Vt0))

	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

	dVtref =  -(Q - Q0)/genData.kQ
	dQ0 = - genData.kReg*genData.ratedPower*(dVt0)
	dP0 = - genData.ratedPower*genData.kReg*(omega)
	return [dElq, dEld, domega, ddelta, dpm, dpPm, dEFD, dvWash, dvPSS, dVtref, dQ0, dP0] #}}}1
