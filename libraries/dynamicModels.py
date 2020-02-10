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
from scipy.optimize import fsolve as fsolve
from copy import deepcopy as copy
import numpy as np
from scipy.integrate import odeint as odeint

# MatplotLib and PyPlot for plotting
import matplotlib.pyplot as pyplot
import matplotlib as mplot
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
mplot.rcParams.update({'font.size': 14})


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
	Eld = 0
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
	Vd = -r*Id - xPd*Iq
	Vt = np.sqrt(Vq**2 + Vd**2)
	
	Tw = genData.Tw
	T1 = genData.T1
	T2 = genData.T2

	EFD = EFD0 + vAVR + vPSS

	dElq = 1/tPdo*(EFD - Elq - (xd - xPd)*Id)
	domega = ( pm - Elq*Iq - (xPd - xPq)*Id*Iq - D*omega )/(2*H)
	ddelta  = omega
	dpm = pPm
	dpPm = (pmSet - pPm*(tG + tT) - pm)/(tG*tT)
	dvAVR = (Ke*(Vt - Vt0) - vAVR)/Te
	dvWash = KPss*domega - vWash/Tw
	dvPSS = (T2*dvWash + vWash - vPSS)/T1

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

def dynamicSimulation(case, disturbanceData, tFinal):
	print(' --> Started dynamic simulation of case \'{}\'...'.format(case.name))
	targetBus, disturbanceAmplitude, disturbanceTime = disturbanceData
	targetBusN = case.getBusNumber(targetBus)

	success = case.runPowerFlow()
	if not success:
		raise Exception(' >>> Dynamic simulation of case {} halted because the initial power flow method did not converge.'.format(case.name))
		
	
	# Building "disturbed case"
	disturbedCase = copy(case)
	disturbedCase.busData[targetBusN].pLoad += np.real(disturbanceAmplitude)
	disturbedCase.busData[targetBusN].qLoad += np.imag(disturbanceAmplitude)

	Yred, C, D, reducedCase = case.reduceMatrixes()
	dYred, dC, dD, disturbedReducedCase = disturbedCase.reduceMatrixes()

	reducedCase.runPowerFlow()

	V = [bus.finalVoltage for bus in case.busData[:case.nGen]]
	I = Yred @ V

	nGen = len(Yred)

	x0 = []
	# Calculating initial generator conditions
	k = 0
	for gen in reducedCase.genData:
		r = gen.ra
		xPd = gen.xPd
		xPq = gen.xPq
		xd = gen.xd
		xq = gen.xq

		V0 = case.busData[gen.busNumber].finalVoltage
		I0 = sum([Yred[gen.busNumber, m]*V[m] for m in range(nGen)])
		
		Vq, Vd = np.real(V0), np.imag(V0)
		Iq, Id = np.real(I0), np.imag(I0)

		gen.vRef = np.sqrt(Vd**2 + Vq**2)

		Eqd = (Vq + 1j*Vd) + (r + 1j*xq)*(Iq + 1j*Id)
		ELq = Vq + r*Iq - xPd*Id

		if gen.modelDepth in [1,2,3,4]: ELd = 0
		else: ELd  = Vd + r*Id - xPq*Iq

		gen.el0 = ELq + 1j*ELd

		delta = np.angle(ELq + 1j*ELd)
		gen.delta0 = delta
		gen.omega0 = 0

		pm0 = ELq*Iq + ELd*Id + (xPd - xPq)*Iq*Id
		if gen.modelDepth == 1: pm0 = ELq*Iq
		gen.pm0 = pm0

		EFD0 = ELq + (xd - xPd)*Id
		gen.efd0 = EFD0
		
		if gen.modelDepth == 1:
			x0.extend([gen.omega0,gen.delta0])
	#		print(SMC([0,delta],Iq + 1j*Id, gen))
		if gen.modelDepth == 2:
			x0.extend([ELq, gen.omega0, gen.delta0])
	#		print(SM1A([ELq,0,delta],Iq + 1j*Id, gen))
		if gen.modelDepth == 3:
			x0.extend([ELq, gen.omega0, gen.delta0, pm0, 0])
		if gen.modelDepth == 4:
			x0.extend([ELq, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0])

	print(x0)
	def dynamicFunction(x,t):
		
		F = []
		targetCase = reducedCase
		Yt = Yred
		if t > disturbanceTime:
			targetCase = disturbedReducedCase
			Yt = dYred

		EL = np.ones(len(Yt), dtype = complex)
		for gen in reducedCase.genData:
			k = gen.busNumber
			if gen.modelDepth == 1:
		#		Iq, Id, Vq, Vd = fsolve(solveI, [1, 0, 1, 0], args = (np.real(gen.el0), np.imag(gen.el0), gen.ra, gen.xPq, gen.xPd, Yt, k, V))
				EL[k] = np.real(gen.el0) + 1j*0

		def solveCurrents(x):
			i = 0
			F = []
			for k in range(nGen):
				gen = reducedCase.genData[k]
				ELq = np.real(EL[k])
				ELd = np.imag(EL[k])

				Iq = x[i]
				Id = x[i+1]
				Vq = ELq - gen.ra*Iq + gen.xPd*Id
				Vd = ELd - gen.ra*Id - gen.xPq*Iq

				F1 = Iq - np.real(sum([Yt[k,m]*V[m] for m in range(nGen)]))
				F2 = Id - np.imag(sum([Yt[k,m]*V[m] for m in range(nGen)]))

				F.extend([F1,F2])
				i += 2
			
			return F

		F0 = [1, 0]*nGen
		sol = fsolve(solveCurrents,F0)
		I = np.ones(nGen, dtype = complex)
		for k in range(nGen):
			I[k] = sol[2*k] + 1j*sol[2*k+1]

		F = []		
		i = 0
		for gen in reducedCase.genData:
			k = gen.busNumber
			if gen.modelDepth == 1:
	#			Iq, Id, Vq, Vd = fsolve(solveI, [1, 0, 1, 0], args = (np.real(gen.el0), np.imag(gen.el0), gen.ra, gen.xPq, gen.xPd, Yt, k, V))
				Iq, Id = np.real(I[k]), np.imag(I[k])
				F.extend(SMC([x[i],x[i+1]], Iq + 1j*Id, gen))
				i += 2
			if gen.modelDepth == 2:
				Iq, Id = np.real(I[k]), np.imag(I[k])
				F.extend(SM1A([x[i], x[i+1], x[i+2]], Iq + 1j*Id, gen))
				i += 3
			if gen.modelDepth == 3:
				Iq, Id = np.real(I[k]), np.imag(I[k])
				F.extend(SM1A_TUR_GOV([x[i], x[i+1], x[i+2], x[i+3], x[i+4]], Iq + 1j*Id, gen))
				i += 5
			if gen.modelDepth == 4:
				Iq, Id = np.real(I[k]), np.imag(I[k])
				F.extend(SM1A_TUR_GOV_AVR_PSS([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7]], Iq + 1j*Id, gen))
				i += 8


				#F.extend(SM1A([x
				#x0.extend([[ELq], [gen.omega0], [gen.delta0]])
		#	if gen.modelDepth == 3:
		#		x0.append([[ELq], [gen.omega0], [gen.delta0], [pm], [0]])
		#	if gen.modelDepth == 4:
		#		x0.append([[ELq], [gen.omega0], [gen.delta0], [pm], [0], [gen.vAVR0], [gen.vWash0], [gen.vPSS0]])
		print(F)
		return F

	tspan = np.linspace(0,tFinal,10001)

	print('HAHA')
	y = odeint(dynamicFunction,x0,tspan)
	fig1, (ax1, ax2) = pyplot.subplots(1,2)

	i = 0
	for gen in reducedCase.genData:
		if gen.modelDepth == 1:
			ax1.plot(tspan,y[:,i], label = reducedCase.busData[gen.busNumber].name)
			ax2.plot(tspan,y[:,i+1], label = reducedCase.busData[gen.busNumber].name)
			i += 2
		if gen.modelDepth == 2:
			ax1.plot(tspan,y[:,i+1], label = reducedCase.busData[gen.busNumber].name)
			ax2.plot(tspan,y[:,i+2], label = reducedCase.busData[gen.busNumber].name)
			i += 3
		if gen.modelDepth == 3:
			ax1.plot(tspan,y[:,i+1], label = reducedCase.busData[gen.busNumber].name)
			ax2.plot(tspan,y[:,i+2], label = reducedCase.busData[gen.busNumber].name)
			i += 5
		if gen.modelDepth == 4:
			ax1.plot(tspan,y[:,i+1], label = reducedCase.busData[gen.busNumber].name)
			ax2.plot(tspan,y[:,i+2], label = reducedCase.busData[gen.busNumber].name)
			i += 8


	ax1.set_ylabel(r'Speed $\omega$')
	ax2.set_ylabel(r'Angle deviation $\delta$')
	ax1.legend()
	ax1.grid(True)
	ax2.legend()
	ax2.grid(True)

	pyplot.show()
