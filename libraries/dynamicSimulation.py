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

import progressbar

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
mplot.rcParams.update({'font.size': 14})

# Definig dynamic simulation routine
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

	# Reducing cases
	Yred, reducedCase = case.reduceMatrixes()
	dYred, disturbedReducedCase = disturbedCase.reduceMatrixes()

	# Calculating initial voltages and currents in phasor form
	V0 = [bus.finalVoltage*np.e**(1j*bus.finalAngle*np.pi/180) for bus in case.busData[:case.nGen]]
	I0 = Yred @ V0

	# deltaQD stores the angles of each machine QD axis in reference to the systems synchronous axis
	deltaQD = [ np.angle(V0[k] + (case.genData[k].ra + 1j*case.genData[k].xq)*I0[k]) for k in range(case.nGen)]

	print(' >>> Voltages calculated at initial time = {}'.format(V0))
	print(' >>> Currents calculated at initial time = {}'.format(I0))
	print(' >>> Axis references calculated at initial time = {}'.format(deltaQD))

	x0 = []
	# Calculating initial generator conditions
	k = 0

	nGen = len(reducedCase.genData)
	EL = np.zeros(case.nGen, dtype = complex)
	for gen in reducedCase.genData:
		#input(' !--> Press any key to continue')
		r = gen.ra
		xPd = gen.xPd
		xPq = gen.xPq

		xd = gen.xd
		xq = gen.xq

		k = gen.busNumber
		print('\n >>> Bus {}'.format(k))
		#print('\n >>> r = {}, xPq =  {}, xPd = {}'.format(r, xPq, xPd))

		VMachine = V0[k]*np.e**(-1j*deltaQD[k])
		IMachine = I0[k]*np.e**(-1j*deltaQD[k])
		Vq, Vd = np.real(VMachine), np.imag(VMachine)
		Iq, Id = np.real(IMachine), np.imag(IMachine)
		print(' >>> Vt = {}'.format(Vq + 1j*Vd))
		#print(' >>> Ib = {}.'.format(Iq + 1j*Id))

		gen.vRef = np.abs(V0[k])
		#print( ' >>> gen VREF = {}, Vt = {}'.format(gen.vRef, np.sqrt(Vq**2 + Vd**2)))
		#print(' >> Set gen VREF = {}'.format(gen.vRef))

		ELq = Vq + r*Iq - xPd*Id
		ELd = Vd + r*Id + xPq*Iq

		EL[k] = ELq + 1j*ELd
		gen.el0 = EL[k]
		
		#print(' >>> EL = {}'.format(EL[k]))

		delta = np.angle(ELq + 1j*ELd)
		gen.delta0 = delta
		gen.omega0 = 0

		pm0 = ELq*Iq + ELd*Id + (xPd - xPq)*Iq*Id
		if gen.modelDepth == 1: pm0 = ELq*Iq
		gen.pm0 = pm0
		#print(' >>> PM0 = {}'.format(pm0))

		EFD0 = ELq + (xd - xPd)*Id
		gen.efd0 = EFD0
		#print(' >>> EFD0 = {}'.format(EFD0))

		gen.vAVR0 = gen.Ke*(np.abs(V0[k]) - gen.vRef)
		dvAVR = (gen.Ke*(np.abs(V0[k]) - gen.vRef) - gen.vAVR0)/gen.Te
		#print(' >>> dvAVR = {}'.format(dvAVR))
		#print(' >>> VRef = {}, V0 = {}'.format(gen.vRef, np.abs(V0[k])))
	
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
			print(' >>> Differential equations for gen at bus \'{}\' at initial time: {}'.format(case.busData[gen.busNumber].name, SM1A_TUR_GOV_AVR_PSS([ELq, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0],Iq + 1j*Id, gen)))

	#print(' >>> Induced voltages at initial point = {}\n'.format(EL))
	#print(' >>> deltaQD at initial point = {}\n'.format(deltaQD))
	#input(' !--> Press any key to continue')

	#print(x0)

	#Imachine = [I0[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	VM = [V0[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	#DIFF = [EL[k] - Vmachine[k] - case.genData[k].ra*Imachine[k] + case.genData[k].xPd*np.imag(Imachine[k]) - 1j*case.genData[k].xPq*np.real(Imachine[k]) for k in range(nGen)]
	#print(DIFF)	

	def solveCurrents(y, *args):
		Y, EL, angleReferences, case = args
		#print(' >>> Passed angleReferences = {}'.format(angleReferences))
		print(' >>> Passed EL = {}'.format(EL))

		i = 0
		nGen  = case.nGen
	
		F = np.zeros(2*nGen, dtype = float)

		V = [y[k] + 1j*y[k + nGen] for k in range(nGen)]

		I = Y @ V

		for k in range(case.nGen):
			# I and V are the current and voltage vectors in the synchronous reference

			ELq = np.real(EL[k])
			ELd = np.imag(EL[k])
		
			Imachine = I[k]*np.e**(-1j*angleReferences[k])
			Vmachine = V[k]*np.e**(-1j*angleReferences[k])

			#print(' >>> Imachine = {}'.format(Imachine))
			Iq = np.real(Imachine)
			Id = np.imag(Imachine)

			Vq = np.real(Vmachine)
			Vd = np.imag(Vmachine)
			
			F[k] = Vq - ELq + gen.ra*Iq - gen.xPd*Id
			F[k + nGen] = Vd - ELd + gen.ra*Id + gen.xPq*Iq
			
		#print(' >>> Calculated F function norm = {}'.format(np.linalg.norm(np.concatenate((I1, I2, V1, V2), axis=0))))
		return F

	print(' >>> solveCurrents initial value = {}'.format(solveCurrents([np.real(x) for x in V0] + [np.imag(x) for x in V0], Yred, EL, deltaQD, case)))

	F0 = [np.real(x) for x in VM] + [np.imag(x) for x in VM]
	sol = fsolve(solveCurrents, F0, args = (Yred, EL, deltaQD, reducedCase))
	print(' >>> solveCurrents initial residual values: {}'.format(solveCurrents(sol, Yred, EL, deltaQD, reducedCase)))

	V = [ sol[k] + 1j*sol[k + nGen] for k in range(nGen)]
	I = Yred @ V
	Vmachine = [V[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	Imachine = [I[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	DIFF = [EL[k] - Vmachine[k] - case.genData[k].ra*Imachine[k] + case.genData[k].xPd*np.imag(Imachine[k]) - 1j*case.genData[k].xPq*np.real(Imachine[k]) for k in range(nGen)]
	print(' >>> solveCurrents initial CALCULATED residual values: {}'.format(DIFF))


	print('\n >>> Initial solveCurrents solution: V = {}\n'.format(V))

	def dynamicFunction(x,t):
		if t > tFinal: pbar.update(tFinal)
		else: pbar.update(t)
		F = []
		targetCase = reducedCase
		Yt = Yred
		if t > disturbanceTime:
			targetCase = disturbedReducedCase
			Yt = dYred

		EL = np.ones(targetCase.nGen, dtype = complex)
		angleReferences = np.ones(targetCase.nGen, dtype = complex)
		i = 0
		for gen in reducedCase.genData:
			k = gen.busNumber
			if gen.modelDepth == 1:
				EL[k] = np.real(elq0) + 1j*0
				angleReferences[k] = deltaQD[k] + x[i+1]
				i += 2
			if gen.modelDepth == 2:
				EL[k] = (x[i] + np.imag(gen.el0))#*np.exp(1j*(deltaQD[k] + x[i+2]))
				angleReferences[k] = deltaQD[k] + x[i+2]
				i += 3
			if gen.modelDepth == 3:
				EL[k] = x[i] + np.imag(gen.el0)
				angleReferences[k] = deltaQD[k] + x[i+2]
				i += 5
			if gen.modelDepth == 4:
				EL[k] = x[i] + np.imag(gen.el0)
				angleReferences[k] = deltaQD[k] + x[i+2]
				i += 8

		#print(' >> Dynamic function EL = {}'.format(EL))
	

		F0 = [1]*(nGen) + [0]*(nGen)
	
		sol = fsolve(solveCurrents,F0, args = (Yt, EL, angleReferences, case))
		V = [ sol[k] + 1j*sol[k + nGen] for k in range(nGen)]
		I = Yt @ V

		I = [ I[k]*np.e**(-1j*angleReferences[k]) for k in range(nGen)]

	#	print(' >>> Induced voltages calculated at time {} = {}'.format(t, EL))	
	#	print(' >>> Current / voltage difference at time {} = {}'.format(t, np.linalg.norm(I - Yt @ V)))
	#	print(' >>> Bus voltages calculated at time {} = {}'.format(t, V))
	#	print(' >>> Bus currents calculated at time {} = {}\n'.format(t,I))

		#print(' >>> Dynamic function I = {}'.format(I))

		F = []		
		i = 0
		for gen in reducedCase.genData:
			k = gen.busNumber
			Imachine = I[k]*np.exp(-1j*deltaQD[k])
			Iq, Id = np.real(I[k]), np.imag(I[k])
			if gen.modelDepth == 1:
				F.extend(SMC([x[i],x[i+1]], Iq + 1j*Id, gen))
				i += 2
			if gen.modelDepth == 2:
				F.extend(SM1A([x[i], x[i+1], x[i+2]], Iq + 1j*Id, gen))
				i += 3
			if gen.modelDepth == 3:
				F.extend(SM1A_TUR_GOV([x[i], x[i+1], x[i+2], x[i+3], x[i+4]], Iq + 1j*Id, gen))
				i += 5
			if gen.modelDepth == 4:
				F.extend(SM1A_TUR_GOV_AVR_PSS([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7]], Iq + 1j*Id, gen))
				i += 8
		#print(F)
		return F

	tspan = np.linspace(0,tFinal,10001)

	print(' --> Integrating differential equations... ')
	pbar = progressbar.ProgressBar(max_value=tFinal)
	y = odeint(dynamicFunction,x0,tspan)
	pbar.finish()
	fig1, (ax1, ax2) = pyplot.subplots(1,2)
	fig2, (ax3, ax4) = pyplot.subplots(1,2)
	fig3, (ax5, ax6) = pyplot.subplots(1,2)

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
			EFD = [gen.efd0 + y[k,i+5] + y[k,i+6] for k in range(len(tspan))] 
			ax3.plot(tspan,y[:,i], label = reducedCase.busData[gen.busNumber].name) 
			ax4.plot(tspan,EFD, label = reducedCase.busData[gen.busNumber].name)
			ax5.plot(tspan, y[:,i+5], label = reducedCase.busData[gen.busNumber].name)
			ax6.plot(tspan, y[:,i+7], label = reducedCase.busData[gen.busNumber].name)
			i += 8

	ax1.set_ylabel(r'Speed $\omega$')
	ax2.set_ylabel(r'Angle deviation $\delta$')
	ax3.set_ylabel(r'Induced voltage $E^\prime_q$')
	ax4.set_ylabel(r'Field voltage $E_{FD}$')
	ax5.set_ylabel(r'AVR voltage $V_{AVR}$')
	ax6.set_ylabel(r'PSS voltage $V_{PSS}$')

	ax1.legend()
	ax1.grid(True)

	ax2.legend()
	ax2.grid(True)

	ax3.legend()
	ax3.grid(True)

	ax4.legend()
	ax4.grid(True)

	ax5.legend()
	ax5.grid(True)

	ax6.legend()
	ax6.grid(True)

	pyplot.show()
