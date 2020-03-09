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

from libraries.dynamicModels import *

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
mplot.rcParams.update({'font.size': 14})

# Definig dynamic simulation routine
def dynamicSimulation(case, disturbanceData, tFinal):
	print(' --> Started dynamic simulation of case \'{}\'...'.format(case.name))
	targetBus, disturbanceAmplitude, disturbanceTime = disturbanceData
	nGen = case.nGen
	nBus = case.nBus
	success = case.runPowerFlow()
	if not success: raise Exception(' >>> Dynamic simulation of case {} halted because the initial power flow method did not converge.'.format(case.name))

	# Reducing case
	Yr, rCase = case.reduceMatrixes()
	targetBusN = rCase.getBusNumber(targetBus)

	# Building Y matrix with new load
	#print('\n\n\n\n\n\n {} \n\n\n\n\n'.format(np.diag([ (np.real(disturbanceAmplitude)  - 1j*np.imag(disturbanceAmplitude))/rCase.busData[targetBusN].finalVoltage**2 if k == targetBusN else 0 for k in range(rCase.nBus)])/rCase.Sb))

	# Coupling coefficients matrix
	print(rCase)

	pinningList = [[0,1,7],[3,2,8],[5,4,6]]

	a = np.zeros((nGen, nGen), dtype = float)
	for k in range(nGen):
		gen = rCase.genData[k]
		for m in range(nGen):
			for i in range(len(pinningList)):
				if k != m and k in pinningList[i] and m in pinningList[i]: a[k,m] = gen.H*np.abs(Yr[k,m])

		a[k,k] = -sum(a[k])

	print(' --> Pinning coefficients matrix:\n{}'.format(a))

	kb = 1e5
	b0 = kb/100	# Generator at bus 0
	b3 = kb/72.5	# Generator at bus 6
	b5 = kb/52.5	# Generator at bus 9
	b = np.array([[b0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,b3,0,0],[0,0,0,0,0,0],[0,0,0,0,0,b5]])

	print(' --> Pinning gains matrix:\n{}'.format(b))

	print(' --> a - b matrix:\n{}'.format(a-b))

	print(' --> EIGENVALUES: {}'.format(np.linalg.eigvals(a-b)))
 

	# Building 'disturbed case'
	dCase = copy(case)
	dCase.name = 'Disturbed ' + dCase.name
	dCase.busData[targetBusN].pLoad += np.real(disturbanceAmplitude)
	dCase.busData[targetBusN].qLoad += np.imag(disturbanceAmplitude)

	dYr, rdCase = dCase.reduceMatrixes()

#	Y = copy(rCase.Y)
#	Y += dYL
#
#	Y1 = Y[0 : nGen, 0 : nGen]
#	Y2 = Y[0 : nGen , nGen : nBus]
#	Y3 = Y[nGen : nBus , 0 : nGen]
#	Y4 = Y[nGen : nBus , nGen : nBus]
#
#	dYr = Y1 - Y2 @ np.linalg.inv(Y4) @ Y3	# dYr is the equivalent reduced matrix of the disturbed case
	#print('\n >>> Disturbed equivalent reduced conductance matrix = \n {}\n'.format(dYr))
	# Calculating initial voltages and currents in phasor form
	V0 = [bus.finalVoltage*np.e**(1j*np.pi/180*bus.finalAngle) for bus in rCase.busData[:rCase.nGen]]
	I0 = Yr @ V0
	#print(' >>> Voltages as obtained through power flow: {}'.format(V0))
	#print(' >>> Currents as obtained through power flow: {}'.format(I0))

	# deltaQD stores the angles of each machine QD axis in reference to the systems synchronous axis
	EQD = [V0[k] + (rCase.genData[k].ra + 1j*rCase.genData[k].xq)*I0[k] for k in range(rCase.nGen)]
	deltaQD = [ np.angle(x) for x in EQD]
	voltageAngles = [np.pi*bus.finalAngle/180 for bus in rCase.busData[:rCase.nGen]]
	#print(' >>> Voltages calculated at initial time = {}'.format(V0))
	#print(' >>> Currents calculated at initial time = {}'.format(I0))
	#print(' >>> Axis references calculated at initial time = {}'.format(deltaQD))

	x0 = []
	# Calculating initial generator conditions
	k = 0

	nGen = len(rCase.genData)
	EL = np.zeros(rCase.nGen, dtype = complex)
	for gen in rCase.genData:
		#input(' !--> Press any key to continue')
		r = gen.ra
		xPd = gen.xPd
		xPq = gen.xPq

		xd = gen.xd
		xq = gen.xq

		k = rCase.getBusNumber(gen.busName)
		print('\n >>> Gen {} at bus \'{}\''.format(k, gen.busName))
		#print('\n >>> r = {}, xPq =  {}, xPd = {}'.format(r, xPq, xPd))

		VMachine = V0[k]*np.e**(-1j*deltaQD[k])
		IMachine = I0[k]*np.e**(-1j*deltaQD[k])
		Vq = np.real(VMachine)
		Vd = np.imag(VMachine)
		Iq = np.real(IMachine)
		Id = np.imag(IMachine)
		print(' >>> Vt = {}'.format(Vq + 1j*Vd))
		print(' >>> Ib = {}.'.format(Iq + 1j*Id))

		gen.vRef = np.abs(V0[k])
		gen.V0 = gen.vRef
		#print( ' >>> gen VREF = {}, Vt = {}'.format(gen.vRef, np.sqrt(Vq**2 + Vd**2)))
		print(' >> Set gen VREF = {}'.format(gen.vRef))

		ELq = Vq + r*Iq - xPd*Id
		ELd = Vd + r*Id + xPq*Iq

		EL[k] = ELq + 1j*ELd
		gen.el0 = EL[k]
		
		print(' >>> EL = {}'.format(EL[k]))

		gen.delta0 = deltaQD[k]
		gen.omega0 = 0
		gen.w0 = gen.omega0

		pm0 = ELq*Iq + ELd*Id + (xPd - xPq)*Iq*Id
		if gen.modelDepth == 1: pm0 = ELq*Iq
		gen.pm0 = pm0
		gen.P0 = pm0

		gen.kP = gen.ratedPower/20

		#print(' >>> PM0 = {}'.format(pm0))

		EFD0 = ELq - (xd - xPd)*Id
		gen.efd0 = EFD0
		#print(' >>> EFD0 = {}'.format(EFD0))

		gen.vAVR0 = -gen.Ke*(np.abs(V0[k]) - gen.vRef)
		#print(' >>> dvAVR = {}'.format(dvAVR))
		#print(' >>> VRef = {}, V0 = {}'.format(gen.vRef, np.abs(V0[k])))

		gen.Q0 = Vd*Iq - Vq*Id
		gen.kQ = 0.1*gen.ratedPower
		gen.kReg = 0.01

		if gen.modelDepth == 1:	x0.extend([gen.omega0,gen.delta0])
		if gen.modelDepth == 2: x0.extend([ELq, gen.omega0, gen.delta0])
		if gen.modelDepth == 3:	x0.extend([ELq, gen.omega0, gen.delta0, pm0, 0])
		if gen.modelDepth in [4,5]: x0.extend([ELq, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0])
		if gen.modelDepth == 6: x0.extend([ELq, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0, gen.vRef])
		if gen.modelDepth in [7,8]: x0.extend([ELq, ELd, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0])
		if gen.modelDepth == 9:	x0.extend([ELq, ELd, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0, gen.vRef, gen.Q0])
		if gen.modelDepth == 10: x0.extend([ELq, ELd, gen.omega0, gen.delta0, pm0, 0])
		if gen.modelDepth == 11: x0.extend([ELq, ELd, gen.omega0, gen.delta0, pm0, 0, gen.vAVR0, gen.vWash0, gen.vPSS0, gen.vRef, gen.Q0, gen.P0])


	#print(' >>> Induced voltages at initial point = {}\n'.format(EL))
	#print(' >>> deltaQD at initial point = {}\n'.format(deltaQD))
	#input(' !--> Press any key to continue')

	#print(x0)

	#Imachine = [I0[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	VM = [V0[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	print(V0)
	print(' >>> Calculated internal induced voltages = {}'.format(EL))
	#DIFF = [EL[k] - Vmachine[k] - case.genData[k].ra*Imachine[k] + case.genData[k].xPd*np.imag(Imachine[k]) - 1j*case.genData[k].xPq*np.real(Imachine[k]) for k in range(nGen)]
	#print(DIFF)	

	print(' >>> solveCurrents initial value = {}'.format(solveCurrents([np.real(x) for x in V0] + [np.imag(x) for x in V0], Yr, EL, deltaQD, rCase)))

	F0 = [np.real(x) for x in VM] + [np.imag(x) for x in VM]
	sol = fsolve(solveCurrents, F0, args = (Yr, EL, deltaQD, rCase))
	print(' >>> solveCurrents initial residual values: {}'.format(solveCurrents(sol, Yr, EL, deltaQD, rCase)))

	V = [ sol[k] + 1j*sol[k + nGen] for k in range(nGen)]
	I = Yr @ V
	Vmachine = [V[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	Imachine = [I[k]*np.e**(-1j*deltaQD[k]) for k in range(nGen)]
	DIFF = [EL[k] - Vmachine[k] - rCase.genData[k].ra*Imachine[k] + rCase.genData[k].xPd*np.imag(Imachine[k]) - 1j*rCase.genData[k].xPq*np.real(Imachine[k]) for k in range(nGen)]
	print(' >>> solveCurrents initial CALCULATED residual values: {}'.format(DIFF))


	print('\n >>> Initial solveCurrents solution: V = {}'.format(V))
	print('\n >>> Powerflow calculated solution V = {}\n'.format(V0))
	def dynamicFunction(x,t, *args):
		Yr, dYr = args
		if t > tFinal: pbar.update(tFinal)
		else: pbar.update(t)
		F = []
		targetCase = rCase
		Yt = Yr
		if t > disturbanceTime:
			Yt = dYr
			targetCase = rdCase

		EL = np.ones(targetCase.nGen, dtype = complex)
		angleReferences = np.zeros(targetCase.nGen)
		i = 0
		for gen in rCase.genData:
			k = rCase.getBusNumber(gen.busName)
			if gen.modelDepth == 1:
				EL[k] = np.real(elq0) + 1j*0
				angleReferences[k] = x[i+1]
				i += 2
			elif gen.modelDepth == 2:
				EL[k] = (x[i] + 1j*np.imag(gen.el0))#*np.exp(1j*(deltaQD[k] + x[i+2]))
				angleReferences[k] = x[i+2]
				i += 3
			elif gen.modelDepth == 3:
				EL[k] = x[i] + 1j*np.imag(gen.el0)
				angleReferences[k] = x[i+2]
				i += 5
			elif gen.modelDepth in [4,5]:
				EL[k] = x[i] + 1j*np.imag(gen.el0)
				angleReferences[k] = x[i+2]
				i += 8
			elif gen.modelDepth == 6:
				EL[k] = x[i] + 1j*np.imag(gen.el0)
				angleReferences[k] = x[i+2]
				i += 9
			elif gen.modelDepth in [7,8]:
				EL[k] = x[i] + 1j*x[i+1]
				angleReferences[k] = x[i+3]
				i += 9
			elif gen.modelDepth == 9:
				EL[k] = x[i] + 1j*x[i+1]
				angleReferences[k] = x[i+3]
				i += 11
			elif gen.modelDepth == 10:
				EL[k] = x[i] + 1j*x[i+1]
				angleReferences[k] = x[i+3]
				i += 6
			elif gen.modelDepth == 11:
				EL[k] = x[i] + 1j*x[i+1]
				angleReferences[k] = x[i+3]
				i += 12

		F0 = [1]*(nGen) + [0]*(nGen)
	
		sol = fsolve(solveCurrents,F0, args = (Yt, EL, angleReferences, targetCase))
		V = [ sol[k] + 1j*sol[k + nGen] for k in range(nGen)]
		I = Yt @ V

		I = [ I[k]*np.e**(-1j*angleReferences[k]) for k in range(nGen)]

		F = []		
		i = 0
		for gen in rCase.genData:
			k = rCase.getBusNumber(gen.busName)
			Iq, Id = np.real(I[k]), np.imag(I[k])
			if gen.modelDepth == 1:
				F.extend(SMC([x[i],x[i+1]], Iq + 1j*Id, gen))
				i += 2
			elif gen.modelDepth == 2:
				F.extend(SM1A([x[i], x[i+1], x[i+2]], Iq + 1j*Id, gen))
				i += 3
			elif gen.modelDepth == 3:
				F.extend(SM1A_TUR_GOV([x[i], x[i+1], x[i+2], x[i+3], x[i+4]], Iq + 1j*Id, gen))
				i += 5
			elif gen.modelDepth == 4:
				F.extend(SM1A_TUR_GOV_AVR_PSS([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7]], Iq + 1j*Id, gen))
				i += 8
			elif gen.modelDepth == 5:
				F.extend(SM1A_TUR_GOV_AVR_PSS_PFQV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7]], Iq + 1j*Id, gen))
				i += 8
			elif gen.modelDepth == 6:
				F.extend(SM1A_TUR_GOV_AVR_PSS_PFQdV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7], x[i+8]], Iq + 1j*Id, gen))
				i += 9
			elif gen.modelDepth == 7:
				F.extend(SM2A_TUR_GOV_AVR_PSS([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7], x[i+8]], Iq + 1j*Id, gen))
				i += 9
			elif gen.modelDepth == 8:
				F.extend(SM2A_TUR_GOV_AVR_PSS_PFQV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7], x[i+8]], Iq + 1j*Id, gen))
				i += 9
			elif gen.modelDepth == 9:
				F.extend(SM2A_TUR_GOV_AVR_PSS_PFQdV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7], x[i+8], x[i+9], x[i+10]], Iq + 1j*Id, gen))
				i += 11
			elif gen.modelDepth == 10:
				F.extend(SM2A_TUR_GOV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5]], Iq + 1j*Id, gen))
				i += 6
			elif gen.modelDepth == 11:
				F.extend(SM2A_TUR_GOV_AVR_PSS_iPFQdV([x[i], x[i+1], x[i+2], x[i+3], x[i+4], x[i+5], x[i+6], x[i+7], x[i+8], x[i+9], x[i+10], x[i+11]], Iq + 1j*Id, gen))
				i += 12

		#print(' >>> Induced voltages calculated at time {} = {}'.format(t, EL))	
		#print(' >>> Current / voltage difference at time {} = {}'.format(t, np.linalg.norm(I - Yt @ V)))
		#print(' >>> Bus voltages calculated at time {} = {}'.format(t, V))
		#print(' >>> Bus currents calculated at time {} = {}\n'.format(t,I))
		#print(' >>> Angle references calculated at time {} = {}\n'.format(t,angleReferences))
		#print(' >>> Dynamic function I = {}'.format(I))
		#print(' >>> Dynamic simulation function: t = {}, F = {}, angleReference = {}'.format(t, F, angleReferences))
		#input()
		return F
	tspan = np.linspace(0,tFinal,10001)

	print(' --> Integrating differential equations... ')
	pbar = progressbar.ProgressBar(max_value=tFinal)
	y = odeint(dynamicFunction,x0,tspan, args = (Yr, dYr))
	pbar.finish()

	# Processing results
	print('\n >>> Simulation done. Processing results...')
	print('  >> Retrieving induced voltages...')

	i = 0
	EL = np.zeros((rCase.nGen, len(tspan)), dtype = complex)
	angleReferences = np.zeros((rCase.nGen, len(tspan)), dtype = float)
	for gen in rCase.genData:
		k = rCase.getBusNumber(gen.busName)
		if gen.modelDepth == 1:
			EL[k] = [np.real(elq0) + 1j*0]*len(tspan)
			angleReferences[k] = y[:,i+1]
			i += 2
		elif gen.modelDepth == 2:
			EL[k] = (y[:,i] + 1j*np.imag(gen.el0))#*np.exp(1j*(deltaQD[k] + x[i+2]))
			angleReferences[k] = y[:,i+2]
			i += 3
		elif gen.modelDepth == 3:
			EL[k] = y[:,i] + 1j*np.imag(gen.el0)
			angleReferences[k] = y[:,i+2]
			i += 5
		elif gen.modelDepth in [4,5]:
			EL[k] = [y[m,i] + 1j*np.imag(gen.el0) for m in range(len(tspan))]
			angleReferences[k] = y[:,i+2]
			i += 8
		elif gen.modelDepth == 6:
			EL[k] = [y[m,i] + 1j*np.imag(gen.el0) for m in range(len(tspan))]
			angleReferences[k] = y[:,i+2]
			i += 9
		elif gen.modelDepth in [7,8]:
			EL[k] = [y[m,i] + 1j*y[m,i+1] for m in range(len(tspan))]
			angleReferences[k] = y[:,i+3]
			i += 9
		elif gen.modelDepth == 9:
			EL[k] = [y[m,i] + 1j*y[m,i+1] for m in range(len(tspan))]
			angleReferences[k] = y[:,i+3]
			i += 11
		elif gen.modelDepth == 10:
			EL[k] = y[:,i] + 1j*y[:,i+1]
			angleReferences[k] = y[:,i+2]
			i += 6
		elif gen.modelDepth == 11:
			EL[k] = [y[m,i] + 1j*y[m,i+1] for m in range(len(tspan))]
			angleReferences[k] = y[:,i+3]
			i += 12

	print('  >> Calculating currents, voltages and output power...')
	V = copy(EL)
	I = copy(EL)
	P = np.zeros((rCase.nGen, len(tspan)), dtype=float)
	Q = np.zeros((rCase.nGen, len(tspan)), dtype=float)
	vRef = np.zeros((rCase.nGen, len(tspan)), dtype=float)
	for i in range(len(tspan)):
		if tspan[i] > disturbanceTime: Yt = dYr
		else: Yt = Yr
		F0 = [1]*(nGen) + [0]*(nGen)
		sol = fsolve(solveCurrents,F0, args = (Yt, EL[:,i], angleReferences[:,i], rCase))
		V[:,i] = [ sol[k] + 1j*sol[k + nGen] for k in range(nGen)]
		I[:,i] = Yt @ V[:,i]

		V[:,i] = [V[k,i]*np.e**(-1j*angleReferences[k,i]) for k in range(nGen)]
		I[:,i] = [I[k,i]*np.e**(-1j*angleReferences[k,i]) for k in range(nGen)]
		P[:,i] = [np.real(V[k,i])*np.real(I[k,i]) + np.imag(V[k,i])*np.imag(I[k,i]) for k in range(nGen)]
		Q[:,i] = [np.imag(V[k,i])*np.real(I[k,i]) - np.real(V[k,i])*np.imag(I[k,i]) for k in range(nGen)]

	i = 0
	print('   >> Calculating AVR reference voltages...')
	for m in range(nGen):
		gen = rCase.genData[m]
		if gen.modelDepth == 1:
			vRef[m,:] = [gen.vRef]*len(tspan)
			i += 2
		elif gen.modelDepth == 2:
			vRef[m,:] = [gen.vRef]*len(tspan)
			i += 3
		elif gen.modelDepth == 3:
			vRef[m,:] = [gen.vRef]*len(tspan)
			i += 5
		elif gen.modelDepth == 4:
			vRef[m,:] = [gen.vRef]*len(tspan)
			i += 8
		elif gen.modelDepth == 5:
			vRef[m,:] = [gen.vRef - (Q[m,t] - gen.Q0)/gen.kQ for t in range(len(tspan))]
			i += 8
		elif gen.modelDepth == 7:
			vRef[m,:] = gen.vRef
			i += 9
		elif gen.modelDepth == 8:
			vRef[m,:] = [gen.vRef - (Q[m,t] - gen.Q0)/gen.kQ for t in range(len(tspan))]
			i += 9
		elif gen.modelDepth == 10:
			vRef[m,:] = [gen.vRef]*len(tspan)
			i += 6


	print('  >> Plotting results...')

	fig1, (ax1,ax2) = pyplot.subplots(1,2)
	fig2, (ax3,ax4) = pyplot.subplots(1,2)
	fig3, (ax5,ax6) = pyplot.subplots(1,2)
	fig4, (ax7,ax8) = pyplot.subplots(1,2)

	i = 0
	for gen in rCase.genData:
		k = rCase.getBusNumber(gen.busName)
		if gen.modelDepth == 1:
			ax1.plot(tspan,y[:,i], label = gen.busName)
			ax2.plot(tspan,y[:,i+1], label = gen.busName)
			i += 2
		if gen.modelDepth == 2:
			ax1.plot(tspan,y[:,i+1], label = gen.busName)
			ax2.plot(tspan,y[:,i+2], label = gen.busName)
			i += 3
		if gen.modelDepth == 3:
			ax1.plot(tspan,y[:,i+1], label = gen.busName)
			ax2.plot(tspan,y[:,i+2], label = gen.busName)
			i += 5
		if gen.modelDepth == 4:
			ax1.plot(tspan,y[:,i+1], label = gen.busName)
			ax2.plot(tspan,y[:,i+2], label = gen.busName)
			EFD = [gen.efd0 + y[k,i+5] + y[k,i+6] for k in range(len(tspan))] 
			ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+5], label = gen.busName)
			ax6.plot(tspan, y[:,i+7], label = gen.busName)
			i += 8
		if gen.modelDepth == 5:
			EFD = [gen.efd0 + y[k,i+5] + y[k,i+7] for k in range(len(tspan))] 
			pmSet = [gen.P0 - gen.kP * y[k, i+1] for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+1], label = gen.busName)

			q = ax2.plot(tspan,pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+3], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan, [np.imag(gen.el0) for t in range(len(tspan))], linestyle = 'dashed', color = q[0].get_color())
			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+5], label = gen.busName)
			ax6.plot(tspan, y[:,i+7], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, vRef[k,:], linestyle ='dashed', label = gen.busName, color = q[0].get_color()) 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, [gen.Q0 for m in range(len(tspan))], linestyle ='dashed', label = gen.busName, color = q[0].get_color()) 

			i += 8
		if gen.modelDepth == 6:
			EFD = [gen.efd0 + y[k,i+5] + y[k,i+7] for k in range(len(tspan))] 
			pmSet = [gen.P0 - gen.kP * y[k, i+1] for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+1], label = gen.busName)

			q = ax2.plot(tspan,pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+3], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan, [np.imag(gen.el0) for t in range(len(tspan))], linestyle = 'dashed', color = q[0].get_color())
			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+5], label = gen.busName)
			ax6.plot(tspan, y[:,i+7], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, y[:,i+8], linestyle ='dashed', label = gen.busName, color = q[0].get_color()) 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, [gen.Q0 for m in range(len(tspan))], linestyle ='dashed', label = gen.busName, color = q[0].get_color()) 

			i += 9

		if gen.modelDepth == 8:
			EFD = [gen.efd0 + y[k,i+6] + y[k,i+8] for k in range(len(tspan))] 
			pmSet = [gen.P0 - gen.kP * y[k, i+2] for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+2], label = gen.busName)

			q = ax2.plot(tspan, pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+4], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan,y[:,i+1], linestyle = 'dashed', color = q[0].get_color(), label = gen.busName) 

			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+6], label = gen.busName)
			ax6.plot(tspan, y[:,i+8], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, vRef[k,:], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference voltage') 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, [gen.Q0 for m in range(len(tspan))], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference reactive power') 

			i += 9
		if gen.modelDepth == 7:
			EFD = [gen.efd0 + y[k,i+6] + y[k,i+8] for k in range(len(tspan))] 
			pmSet = [gen.P0 for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+2], label = gen.busName)

			q = ax2.plot(tspan, pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+4], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan,y[:,i+1], linestyle = 'dashed', color = q[0].get_color(), label = gen.busName) 

			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+6], label = gen.busName)
			ax6.plot(tspan, y[:,i+8], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, vRef[k,:], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference voltage') 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, [gen.Q0 for m in range(len(tspan))], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference reactive power') 

			i += 9

		if gen.modelDepth == 9:
			EFD = [gen.efd0 + y[k,i+6] + y[k,i+8] for k in range(len(tspan))] 
			pmSet = [gen.P0 - gen.kP * y[k, i+2] for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+2], label = gen.busName)

			q = ax2.plot(tspan, pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+4], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan,y[:,i+1], linestyle = 'dashed', color = q[0].get_color(), label = gen.busName) 

			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+6], label = gen.busName)
			ax6.plot(tspan, y[:,i+8], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, y[:,i+9], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference voltage') 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, y[:,i+10], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference reactive power') 

			i += 11
		if gen.modelDepth == 10:
			EFD = [gen.efd0 for k in range(len(tspan))] 
			pmSet = [gen.pm0 for k in range(len(tspan))]
			ax1.plot(tspan,y[:,i+2], label = gen.busName)

			q = ax2.plot(tspan, pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+4], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan,y[:,i+1], linestyle = 'dashed', color = q[0].get_color(), label = gen.busName) 

			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, [0 for k in range(len(tspan))], label = gen.busName)
			ax6.plot(tspan, [0 for k in range(len(tspan))], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, [gen.vRef for k in range(len(tspan))], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference voltage') 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, [gen.Q0 for m in range(len(tspan))], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference reactive power') 

			i += 6
		if gen.modelDepth == 11:
			EFD = [gen.efd0 + y[k,i+6] + y[k,i+8] for k in range(len(tspan))] 
			pmSet = y[:,i+11] - gen.kP * y[:,i+2]
			ax1.plot(tspan,y[:,i+2], label = gen.busName)

			q = ax2.plot(tspan, pmSet, linestyle = 'dashed', label = gen.busName)
			ax2.plot(tspan, y[:, i+4], label = gen.busName, color = q[0].get_color())

			q = ax3.plot(tspan,y[:,i], label = gen.busName) 
			ax3.plot(tspan,y[:,i+1], linestyle = 'dashed', color = q[0].get_color(), label = gen.busName) 

			ax4.plot(tspan,EFD, label = gen.busName)
			ax5.plot(tspan, y[:,i+6], label = gen.busName)
			ax6.plot(tspan, y[:,i+8], label = gen.busName)
			
			q = ax7.plot(tspan, [np.abs(V[k,m]) for m in range(len(tspan))], label = gen.busName) 
			ax7.plot(tspan, y[:,i+9], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference voltage') 

			q = ax8.plot(tspan, Q[k,:], label = gen.busName) 
			ax8.plot(tspan, y[:,i+10], linestyle ='dashed', color = q[0].get_color(), label = gen.busName + ' reference reactive power') 

			i += 12

	ax1.set_ylabel(r'Speed $\omega$')
	ax2.set_ylabel(r'Mechanical power $P_m$')
	ax3.set_ylabel(r'Induced voltage $E^\prime$')
	ax4.set_ylabel(r'Field voltage $E_{FD}$')
	ax5.set_ylabel(r'AVR voltage $V_{AVR}$')
	ax6.set_ylabel(r'PSS voltage $V_{PSS}$')
	ax7.set_ylabel(r'Reference voltage')
	ax8.set_ylabel(r'Output reactive power $Q$')

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

	ax7.legend()
	ax7.grid(True)

	ax8.legend()
	ax8.grid(True)

	print(' >>> Results processed. Plotting...\n')

	pyplot.show()
