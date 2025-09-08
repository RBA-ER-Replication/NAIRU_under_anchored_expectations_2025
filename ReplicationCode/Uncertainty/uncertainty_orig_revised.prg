'Estimate mean-squared error (MSE) of NAIRU estimates, accounting for both filtering and parameter uncertainty.
'Matthew Read, Economic Research, December 2023
' Extended to probability and modularity by Ballantyne Jan 2024


'Specify name of state-space model object
%modname = "model_"+varname+"_huber"  ' Use Huber-White covariance version of model
%v = slack_measure
%varname = varname

smpl {est_start} {est_end}

' Priors
if AR_inf_error_orig = 1 then
vector(2) snairu
snairu.fill 1.5, 0.5						   	
sym(2) vnairu
vnairu.fill 3.5
else
vector(1) snairu
snairu.fill 1.5						   	
sym(1) vnairu
vnairu.fill 3.5		
endif

'Save MLE of coefficients in vector
vector coefs_mle = {%modname}.@coefs
scalar nparams = coefs_mle.@rows

'Re-estimate model via maximum likelihood but with robust estimate of variance-covariance matrix
'{%modname}.ml(m=500, cov=huber, showopts)

'Save variance-covariance matrix of MLE in symbolic matrix
sym cov_mle = {%modname}.@coefcov
sym cov_ic = @cholesky(cov_mle)

'Draw vectors from a multivariate normal distribution centered at MLE with variance-covariance matrix equal to estimated variance-covariance matrix of coefficients
matrix(K, coefs_mle.@rows) draws 'Storage
for !i=1 to K

	vector tmp = coefs_mle + @rmvnorm(cov_mle) 'Draw of coefficients
	rowplace(draws, @transpose(tmp), !i) 'Save in matrix

next
delete tmp

' Extract list of coefficients in the correct order of covar matrix
freeze(tmp) {%modname}.coefcov
for !i=1 to nparams
	string coefname_{!i} = tmp(!i+2,1)
next
delete tmp

'For each draw of coefficients, compute fixed-interval smoother estimates of state(s) and their variance-covariance matrix.
matrix(@obssmpl, K) draws_{%varname}_sm
matrix(@obssmpl,K) draws_{%varname}_smse
' Extend if making partial quarter probability 
if partial = "yes" then
	matrix(@obssmpl+1, K) probs_{%varname}
else 
	matrix(@obssmpl, K) probs_{%varname}
endif

' Copy model 
copy %modname model_draw

for !i=1 to K

	'Replace existing coefficients (from estimated model) with random draws obtained above
	for !j=1 to nparams
		{coefname_{!j}} = draws(!i,!j)
	next

	'Use maximum likelihood command (since for whatever reason Eviews won't let me generate the states otherwise), but with max number of iterations (m) set to zero. This means that the smoother will be calculated at the values saved in beta, delta, etc.
	model_draw.ml(m=0,showopts)

	'Generate fixed-interval smoother estimate of state(s) and standard errors (i.e. square root of filter variance) and save to matrix
	model_draw.makestates(t=smooth) *_sm
	stomna(nairu_sm, tmp)
	colplace(draws_{%varname}_sm, tmp, !i)
	model_draw.makestates(t=smoothse) *_smse
	stomna(nairu_smse, tmp)
	colplace(draws_{%varname}_smse, tmp, !i)

	' Make probability of unemp > NAIRU

	' Extend if making partial quarter probability 
	if partial = "yes" then
		smpl {est_start} {est_end}
		series nairu_sm_ext = nairu_sm
		smpl {est_end}+1 {est_end}+1
		series nairu_sm_ext = nairu_sm(-1)  ' extend nairu
		series nairu_smse =  nairu_smse(-1) ' extend rmse
		smpl {est_start} {est_end}+1 ' extended date carries through
		' Make distance from NAIRU
		series dist = {%v}-nairu_sm_ext
	else 
		' Make distance from NAIRU
		series dist = {%v}-nairu_sm
	endif
	
	' Make uncertainty interval and calculate probability based on normal
	series stand = dist/nairu_smse
	series prob = @cnorm(stand) ' The probability that u/util is above NAIRU of model, effectively like a one-sided z test

	' Fill matrix
	stomna(prob, tmp)
	colplace(probs_{%varname}, tmp, !i)

	' Ensure sample reset
	smpl {est_start} {est_end}


next

'Compute average of variance of state estimate across draws (at each time period)
'Square standard errors/RMSE to get smoother variance at each draw
matrix draws_smvar = @emult(draws_{%varname}_smse, draws_{%varname}_smse)
'Average smoother variance across draws
vector smoother_var_mat = @cmean(@transpose(draws_smvar))

'Compute squared deviations between smoothed estimates at each draw and at MLE of parameters
matrix(@obssmpl, K) sm_orig
for !i=1 to K
	stomna(nairu_{%varname}_sm, tmp)
	colplace(sm_orig, tmp, !i)
next

matrix sm_diff = draws_{%varname}_sm - sm_orig
matrix sm_diff_sq = @emult(sm_diff, sm_diff)
vector par_var_mat = @cmean(@transpose(sm_diff_sq))

mtos(smoother_var_mat, smoother_var_{%varname})
mtos(par_var_mat, par_var_{%varname})

' Within sampler probability average
if partial = "yes" then
	smpl {est_start} {est_end}+1
endif
vector probs_avg_{%varname}_vec = @cmean(@transpose(probs_{%varname}))
mtos(probs_avg_{%varname}_vec, probs_avg_{%varname})

