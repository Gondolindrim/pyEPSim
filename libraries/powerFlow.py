
def powerFlow(case,absTol,deltaMax,maxIter)
	print(' --> Beggining power flow...')
	# Initial guess: flat start
	V = .9*np.ones(case.nBus)
	theta = np.zeros(case.nBus)
	P = np.diag(array([br.pGen - br.pLoad for br in case.branchData]))
	Q = np.diag(array([br.qGen - br.qLoad for br in case.branchData]))
	isP = np.eye(case.nBus)
	isV = np.zeros(case.nBus)

	numP = int(isP.sum())
	numV = int(isV.sum())


	itCount = 0	



	# (5.6) Start time counte6
	tstart = time.time()

	# -------------------------------------------------
	# (10) STARTING NUMERICAL METHOD FOR POWER FLOW {{{1
	# -------------------------------------------------

	while(True):
	# (6) STARTING ITERATIONS

		# (6.1) Increasing iteration counter
		itCount += 1
		
		# (6.2) Printing iteration report
		if verbose > 0: print('\n ==> Iteration #{0:3.0f} '.format(itCount) + '-'*50)
		
		# (6.3) Calculating mF.Jacobian
		H = mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV)
		H = np.delete(H,0,axis=0) # Removing slack bar angle derivatives for it is known

		# (6.5) Calculating state update
		deltaSLC = np.delete(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV),0,axis=0)
		dX = inv(H) @ deltaSLC

		# (6.6) Updating V and theta
		theta[1:] += dX[0:nBus-1,0]
		V += dX[nBus-1:,0]
		
		# (6.8) Printing iteration results
		if verbose > 0:
			print(' --> |dX| = {0}\n --> dV = {1}\n --> dTheta = {2}'.format(norm(dX),transpose(dX[nBus-1:,0]),dX[0:nBus-1,0]))
		if verbose > 1:
			print('\n --> J = \n{0},\n\n --> r = Z - h = \n{2}*\n{1}'.format(	mF.Jac(V,theta,K,a,y,Y,bsh,isP,isV),
												(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))/norm(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV)),
												norm(mF.Z(P,Q,isP,V,isV) - mF.h(V,theta,K,a,y,Y,bsh,isP,isV))))

		# (6.9) Testing for iteration sucess or divergence
		if norm(dX) < absTol: # (6.8.1) If success
			print('\n --> Sucess!')
			print('\n --> Numerical method stopped at iteration {1} with |deltaX| = {0}.'.format(norm(dX),itCount))
			break	

		if norm(dX) > deltaMax: #(6.8.2) If diverted
			print('\n --> Error: the solution appears to have diverged on iteration number {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break
		
		if itCount > maxIter - 1: # (6.8.3) If reached maximum iteration
			print(' --> Error: the solution appears to have taken too long. Final iteration: {0}: |deltaX| = {1}.'.format(itCount,norm(dX)))
			break

		# (6.10) Pausing for each iteration
		if verbose>1:
			pause('\n --> Program paused for next iteration. ')


	# (6.11) Calculating elapsed time
	elapsed = time.time() - tstart

	# (6.13) Printing convergence results
	print(' --> Final result:\n\n	theta = {0} \n\n	V     = {1}'.format(theta,V))

	if verbose > 1: print('\n --> Residual = \n{1}\n\n --> Elapsed time: {0} s'.format(elapsed,r))
