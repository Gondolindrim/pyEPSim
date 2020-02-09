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
def SM2(x,I,genData):
	#F = np.zeros(2)
	omega = x[0]
	delta = x[1]

	H = genData[3]
	D = genData[4]
	pm = genData[17]	# Mechanical power is constant equal to its initial value
	Elq0 = genData[19]
	Iq = np.real(I)

	omegaP = ( pm - Elq0*Iq - D*omega )/(2*H)	# omega = x[k+1]
	deltaP = omega
	return [omegaP, deltaP]

# Synchronous machine one-axis (third-order) model
def SM1A(x,I,genData):
	EPq = x[0]
	omega = x[1]
	delta = x[2]

	H = genData[3]
	D = genData[4]
	EFD0 = genData[18]
	tPdo = genData[10]
	xPd = genData[8]
	xPq = genData[12]
	xd = genData[7]
	H = genData[2]
	D = genData[3]
	pm = genData[17]	# Mechanical power is constant equal to its initial value
	Iq = np.real(I)
	Id = np.imag(I)

	dElq = 1/tPdo*(EFD0 - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta = omega

	return [dElq, domega, ddelta]

# Synchronous machine one-axis (third-order) model with turbine and governor equations
def SM1A_TUR_GOV(x,I,genData):
	Elq = x[0]
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	
	H = genData[3]
	D = genData[4]
	pm* = genData[17]	# Initial value of Pm setpoint is the initial power
	EFD0 = genData[18]
	tPdo = genData[10]
	xPd = genData[8]
	xPq = genData[12]
	xd = genData[7]
	tG = genData[30]
	tT = genData[31]
	Iq = np.real(I)
	Id = np.imag(I)

	dElq = 1/tPdo*(EFD0 - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpm* = (pm* - pPm*(tG + tT) - pm)/(tG*tT)

	return [dElQ, domega, ddelta, dpm, dpm*]

# Synchronous machine one-axis (third-order) model with turbine, governor, AVR and PSS equations
def SM1A_TUR_GOV_AVR_PSS(x,I,genData):
	Elq = x[0]
	omega = x[1]
	delta = x[2]
	pm = x[3]
	pPm = x[4]
	vAVR = x[5]
	vWash = x[6]
	vPSS = x[7]
	
	H = genData[3]
	D = genData[4]
	r = genData[5]
	pm* = genData[17]	# Initial value of Pm setpoint is the initial power
	EFD0 = genData[18]
	tPdo = genData[10]
	xPd = genData[8]
	xPq = genData[12]
	xd = genData[7]
	tG = genData[30]
	tT = genData[31]
	Vt0 = genData[34]
	Iq = np.real(I)
	Id = np.imag(I)
	Vtq = Elq - r*Iq + xPd*Id
	Vtd = -rId - xPq*Iq
	Vt = np.sqrt(Vtq**2 + Vtd**2)
	Tw = genData[27]
	T1 = genData[28]
	T2 = genData[29]

	dElq = 1/tPdo*(EFD0 - Elq + (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpm* = (pm* - pPm*(tG + tT) - pm)/(tG*tT)
	dvAVR = (Ke*(Vt - Vt0) - vAVR)/Te
	dvWash = KPss/(2*H)*( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega ) - vWash/Tw
	dvPSS = T2/T1*( KPss/(2*H)*( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega ) +  vWash*(1/T2 - 1/Tw) - 1/T2*vPSS)
	return [dElQ, domega, ddelta, dpm, dpm*]
