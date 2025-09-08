'===========================================================
' AER Policy Forum 2025: The NAIRU under anchored inflation expectations
' Model program - original NAIRU model with new breaks and COVID adjustments

' Ballantyne & Cusbert 2025
' Leveraging Tom Cusbert (https://www.rba.gov.au/publications/bulletin/2017/jun/2.html), and further developed by Mike Read, and later Alex Ballantyne. 
'===========================================================
 
' Include sspace check
include check_sspace.prg

'--------------------------------------------------------------------------------
' 1. Model setup
'--------------------------------------------------------------------------------


'toggle for AR error term in inflation equation
scalar AR_inf_error_orig = 1

' Loop of unemp and util
for %v u' util

' Model name
if %v == "util" then
	%mname = "orig_util"
else
if  %v == "under" then
	%mname = "orig_under"
else
	%mname = "orig"
endif
endif


' Coefficient vectors
coef(100) beta = 0.1 'use for lags of prices and wages for LHS prices eqs
coef(100) phi = 0.1 'use for lags of prices and wages for LHS prices eqs
coef(100) gamma = 0.1 'gamma(1) is inflation expectations parameter for prices, gamma(2) is for wages
coef(100) delta = 0.1 'delta(2) is inflation expectations parameter for prices, gamma(2) is for wages
coef(100) omega = 0.1 'for speed limit term
coef(100) psi = 0.1 'for import prices
coef(200) sigma =0.5 'for error variances
coef(100) xi = 0 'for constant correlations
coef(100) zeta = 0 'for COVID adjustments applied via toggle
c=0


' Priors
if AR_inf_error_orig = 1 then
	vector(2) sprior_{%mname}
	sprior_{%mname}.fill 1.5, 0.5						   	
	sym(2) vprior_{%mname}
	for !i=1 to 2
		vprior_{%mname}(!i,!i) = 3.5
	next
else
	vector(1) sprior_{%mname}
	sprior_{%mname}.fill 1.5						   	
	sym(1) vprior_{%mname}
	vprior_{%mname}.fill 3.5		
endif
'--------------------------------------------------------------------------------
'2. Estimate model
'--------------------------------------------------------------------------------

'Set sample
smpl {est_start} {est_end}

'2.1 OLS to get parameter starting vals

gamma=-.5

equation pc_prices.ls infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) +gamma(1)*({%v}-c(1))/{%v}+(1-IT)*omega(1)*d({%v}(-1))/{%v}+psi(1)*(pmye(-1) - infx(-1))+c(2)*oil_dummy*d(log(oil(-2)))

equation pc_ulc.ls ulcq= delta(3)*infxq+phi(1)*infq(-1)+(1-phi(1)-delta(3))*infq(-2)+gamma(2)*({%v}-c(1))/{%v}+omega(2)*d({%v}(-1))/{%v}+c(3)*oil_dummy*d(log(oil(-2)))

'2.2 State space model
sspace model_{%mname}

if AR_inf_error_orig =1 then
model_{%mname}.append @signal infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) + gamma(1)*({%v}-nairu)/{%v} + (1-IT)*omega(1)*d({%v}(-1))/{%v} + psi(1)*(pmye(-1) - infx(-1)) + c(2)*oil_dummy*d(log(oil(-2))) + AR_E
else
model_{%mname}.append @signal infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) + gamma(1)*({%v}-nairu)/{%v} + (1-IT)*omega(1)*d({%v}(-1))/{%v} + psi(1)*(pmye(-1) - infx(-1)) + c(2)*oil_dummy*d(log(oil(-2))) + [ename = e1, VAR = (SIGMA(1)*(1-IT)+IT*SIGMA(2)+COVID*zeta(20))^2] 
endif

' New spec, IT break and only one lag of infq 
model_{%mname}.append @signal  ulcq = IT*infxq +  (1-IT)*delta(4)*infxq + (1-IT)*(1-delta(4))*infq(-1) + gamma(2)*({%v}-nairu)/{%v} + omega(2)*d({%v}(-1))/{%v} + c(3)*oil_dummy*d(log(oil(-2))) + [ename = e2, var = (sigma(3)*(1-NAV)+NAV*sigma(4))^2] 


model_{%mname}.append @state nairu = nairu(-1) + [ename = e3, VAR = (SIGMA(10)*(1-IT)+IT*SIGMA(11))^2]

model_{%mname}.append @mprior sprior_{%mname}
model_{%mname}.append @vprior vprior_{%mname}

if AR_inf_error_orig=1 then
model_{%mname}.append @state AR_E = COVID_inf_error*zeta(99)*AR_E(-1) + [ename = e1, VAR = (SIGMA(1)*(1-IT)+IT*SIGMA(2))^2] 
endif


model_{%mname}.ml(m=500, showopts)  ' Don't use Huber-White for baseline estimation as it allows estimation of misspecified model (e.g. multicolinearity)
' Make Huber-White covariance version of model at basic MLE point estimates
copy model_{%mname} model_{%mname}_huber
model_{%mname}_huber.ml(m=0, cov=huber, showopts)
' Check convergence - will throw a dialogue box if found
call check_sspace("model_"+%mname)

' Save coefficients for use in uncertainty code
vector coef_{%mname}_beta = beta
vector coef_{%mname}_c = c
vector coef_{%mname}_delta = delta
vector coef_{%mname}_gamma = gamma
vector coef_{%mname}_omega = omega
vector coef_{%mname}_phi = phi
vector coef_{%mname}_psi = psi
vector coef_{%mname}_sigma = sigma
vector coef_{%mname}_xi = xi
vector coef_{%mname}_zeta = zeta

' Create smoothed nairu and standard errors 
model_{%mname}.makestates(t=smooth) *_{%mname}_sm
model_{%mname}.makestates(t=smoothse) *_{%mname}_smse
model_{%mname}.makestates(t=filt) *_{%mname}_filt

' End loop over unemp/util
next


