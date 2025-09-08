'===========================================================
' AER Policy Forum 2025: The NAIRU under anchored inflation expectations
' Model program - original NAIRU model with new breaks, precovid sample

' Ballantyne & Cusbert 2025
' Leveraging Tom Cusbert (https://www.rba.gov.au/publications/bulletin/2017/jun/2.html), and further developed by Mike Read, and later Alex Ballantyne. 
'===========================================================
 

' Include sspace check
include check_sspace.prg

'--------------------------------------------------------------------------------
' 1. Model setup
'--------------------------------------------------------------------------------


' Loop of unemp and util
for %v u 

' Model name
if %v == "util" then
	%mname = "orig_precovid_util"
else
if  %v == "under" then
	%mname = "orig_precovid_under"
else
	%mname = "orig_precovid"
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
	vector(1) sprior_{%mname}
	sprior_{%mname}.fill 1.5						   	
	sym(1) vprior_{%mname}
	vprior_{%mname}.fill 3.5		

'--------------------------------------------------------------------------------
'2. Estimate model
'--------------------------------------------------------------------------------

'Set sample
smpl {est_start} 2019Q4

'2.1 OLS to get parameter starting vals

gamma=-.5

equation pc_prices.ls infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) +gamma(1)*({%v}-c(1))/{%v}+(1-IT)*omega(1)*d({%v}(-1))/{%v}+psi(1)*(pmye(-1) - infx(-1))+c(2)*oil_dummy*d(log(oil(-2)))

equation pc_ulc.ls ulcq= delta(3)*infxq+phi(1)*infq(-1)+(1-phi(1)-delta(3))*infq(-2)+gamma(2)*({%v}-c(1))/{%v}+omega(2)*d({%v}(-1))/{%v}+c(3)*oil_dummy*d(log(oil(-2)))

'2.2 State space model
sspace model_{%mname}

model_{%mname}.append @signal infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) + gamma(1)*({%v}-nairu)/{%v} + (1-IT)*omega(1)*d({%v}(-1))/{%v} + psi(1)*(pmye(-1) - infx(-1)) + c(2)*oil_dummy*d(log(oil(-2))) + [ename = e1, VAR = (SIGMA(1)*(1-IT)+IT*SIGMA(2))^2] 

model_{%mname}.append @signal  ulcq = IT*infxq +  (1-IT)*delta(4)*infxq + (1-IT)*(1-delta(4))*infq(-1) + gamma(2)*({%v}-nairu)/{%v} + omega(2)*d({%v}(-1))/{%v} + c(3)*oil_dummy*d(log(oil(-2))) + [ename = e2, var = (sigma(3)*(1-NAV)+NAV*sigma(4))^2] 

model_{%mname}.append @state nairu = nairu(-1) + [ename = e3, VAR = (SIGMA(10)*(1-IT)+IT*SIGMA(11))^2]

model_{%mname}.append @mprior sprior_{%mname}
model_{%mname}.append @vprior vprior_{%mname}


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



'--------------------------------------------------------------------------------
' B1. Model setup - UNRESTRICTED VERSION
'--------------------------------------------------------------------------------


' Loop of unemp and util
for %v u 'util

' Model name
if %v == "util" then
	%mname = "orig_precovid_unres_util"
else
if  %v == "under" then
	%mname = "orig_precovid_unres_under"
else
	%mname = "orig_precovid_unres"
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
	vector(1) sprior_{%mname}
	sprior_{%mname}.fill 1.5						   	
	sym(1) vprior_{%mname}
	vprior_{%mname}.fill 3.5		

'--------------------------------------------------------------------------------
'B2. Estimate model
'--------------------------------------------------------------------------------

'Set sample
smpl {est_start} 2019Q4

'2.1 OLS to get parameter starting vals

gamma=-.5

equation pc_prices.ls infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) +gamma(1)*({%v}-c(1))/{%v}+(1-IT)*omega(1)*d({%v}(-1))/{%v}+psi(1)*(pmye(-1) - infx(-1))+c(2)*oil_dummy*d(log(oil(-2)))

equation pc_ulc.ls ulcq= delta(3)*infxq+phi(1)*infq(-1)+(1-phi(1)-delta(3))*infq(-2)+gamma(2)*({%v}-c(1))/{%v}+omega(2)*d({%v}(-1))/{%v}+c(3)*oil_dummy*d(log(oil(-2)))

'2.2 State space model
sspace model_{%mname}

model_{%mname}.append @signal infq= IT*delta(1)*infxq+(1-IT)*delta(2)*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) + gamma(1)*({%v}-nairu)/{%v} + (1-IT)*omega(1)*d({%v}(-1))/{%v} + psi(1)*(pmye(-1) - infx(-1)) + c(2)*oil_dummy*d(log(oil(-2))) + [ename = e1, VAR = (SIGMA(1)*(1-IT)+IT*SIGMA(2))^2] 

model_{%mname}.append @signal  ulcq = IT*delta(3)*infxq + (1-IT)*delta(4)*infxq + (1-IT)*phi(4)*infq(-1) + gamma(2)*({%v}-nairu)/{%v} + omega(2)*d({%v}(-1))/{%v} + c(3)*oil_dummy*d(log(oil(-2))) + [ename = e2, var = (sigma(3)*(1-NAV)+NAV*sigma(4))^2] 

model_{%mname}.append @state nairu = nairu(-1) + [ename = e3, VAR = (SIGMA(10)*(1-IT)+IT*SIGMA(11))^2]

model_{%mname}.append @mprior sprior_{%mname}
model_{%mname}.append @vprior vprior_{%mname}


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





'Gather break evidence together


freeze(wald_break_nairu_var) model_orig_precovid_huber.wald @abs(sigma(10))=@abs(sigma(11))
freeze(wald_break_inf_var) model_orig_precovid_huber.wald @abs(sigma(1))=@abs(sigma(2))
freeze(wald_break_ulc_var) model_orig_precovid_huber.wald @abs(sigma(3))=@abs(sigma(4))

freeze(wald_break_inf_infx) model_orig_precovid_huber.wald (1-BETA(1)-BETA(2))=(1-BETA(1)-BETA(2)-BETA(3)-BETA(6))
freeze(wald_break_ulc_infx) model_orig_precovid_huber.wald delta(4)=1


'Breaks in the Phillips curve? Estimate model with extra break for PC slope at IT
' This version also holds the lag structure constant in pre- and post-IT for comparability, and leaves speed limit in.

' Coefficient vectors
coef(100) beta = 0.1 'use for lags of prices and wages for LHS prices eqs
coef(100) phi = 0.1 'use for lags of prices and wages for LHS prices eqs
coef(100) gamma = -0.4 'gamma(1) is inflation expectations parameter for prices, gamma(2) is for wages
coef(100) delta = 0.1 'delta(2) is inflation expectations parameter for prices, gamma(2) is for wages
coef(100) omega = 0.1 'for speed limit term
coef(100) psi = 0.1 'for import prices
coef(200) sigma =0.4 'for error variances
coef(100) xi = 0 'for constant correlations
coef(100) zeta = 0 'for COVID adjustments applied via toggle
c=0
	%mname = "orig_precovid_pcbreak"
%v ="u" 

' Priors
	vector(1) sprior_{%mname}
	sprior_{%mname}.fill 1.5						   	
	sym(1) vprior_{%mname}
	vprior_{%mname}.fill 3.5		


gamma=-.5

equation pc_prices.ls infq= IT*(1-beta(1)-beta(2))*infxq+(1-IT)*(1-beta(1)-beta(2)-beta(3)-beta(6))*infxq + beta(1)*infq(-1) + beta(2)*infq(-2) + (1-IT)*beta(3)*infq(-3) + (1-IT)*beta(6)*(ulcq(-1)) +gamma(1)*({%v}-c(1))/{%v}+(1-IT)*omega(1)*d({%v}(-1))/{%v}+psi(1)*(pmye(-1) - infx(-1))+c(2)*oil_dummy*d(log(oil(-2)))

equation pc_ulc.ls ulcq= delta(3)*infxq+phi(1)*infq(-1)+(1-phi(1)-delta(3))*infq(-2)+gamma(2)*({%v}-c(1))/{%v}+omega(2)*d({%v}(-1))/{%v}+c(3)*oil_dummy*d(log(oil(-2)))



'2.2 State space model
sspace model_{%mname}

model_{%mname}.append @signal infq= (IT)*BETA(11)*INFQ(-1) + (IT)*BETA(21)*INFQ(-2)+(IT)*BETA(31)*INFQ(-3) + (IT)*BETA(61)*(ULCQ(-1))+ IT*(1-BETA(11)-BETA(21)-BETA(31)-BETA(61))*INFXQ+(1-IT)*(1-BETA(1)-BETA(2)-BETA(3)-BETA(6))*INFXQ + (1-IT)*BETA(1)*INFQ(-1) + (1-IT)*BETA(2)*INFQ(-2) + (1-IT)*BETA(3)*INFQ(-3) + (1-IT)*BETA(6)*(ULCQ(-1)) +   (1-it)*GAMMA(1)*({%v}-NAIRU)/{%v}+IT*GAMMA(11)*({%v}-NAIRU)/{%v} + (IT)*omega(1)*d({%v}(-1))/{%v} + (1-IT)*omega(11)*d({%v}(-1))/{%v} + psi(1)*(pmye(-1) - infx(-1)) + c(2)*oil_dummy*d(log(oil(-2))) + [ename = e1, VAR = (SIGMA(1)*(1-IT)+IT*SIGMA(2))^2] 




model_{%mname}.append @signal  ulcq = IT*delta(3)*infxq + (IT)*(1-delta(3))*infq(-1)  +  (1-IT)*delta(4)*infxq + (1-IT)*(1-delta(4))*infq(-1) +  (1-it)*GAMMA(2)*({%v}-NAIRU)/{%v}+IT*GAMMA(21)*({%v}-NAIRU)/{%v}  + (1-IT)*omega(2)*d({%v}(-1))/{%v}+(IT)*omega(22)*d({%v}(-1))/{%v} + c(3)*oil_dummy*d(log(oil(-2))) + [ename = e2, var = (sigma(3)*(1-NAV)+NAV*sigma(4))^2] 

model_{%mname}.append @state nairu = nairu(-1) + [ename = e3, VAR = (SIGMA(10)*(1-IT)+IT*SIGMA(11))^2]

model_{%mname}.append @mprior sprior_{%mname}
model_{%mname}.append @vprior vprior_{%mname}


model_{%mname}.ml(m=500, showopts)  ' Don't use Huber-White for baseline estimation as it allows estimation of misspecified model (e.g. multicolinearity)
' Make Huber-White covariance version of model at basic MLE point estimates
copy model_{%mname} model_{%mname}_huber
model_{%mname}_huber.ml(m=0, cov=huber, showopts)

' Check convergence - will throw a dialogue box if found

call check_sspace("model_"+%mname)


freeze(wald_allbreak_inf_infx) model_{%mname}_huber.wald (1-BETA(11)-BETA(21)-BETA(31)-BETA(61))=(1-BETA(1)-BETA(2)-BETA(3)-BETA(6))
freeze(wald_allbreak_ulc_infx) model_{%mname}_huber.wald delta(4)=delta(3)

freeze(wald_allbreak_inf_pc)  model_orig_precovid_pcbreak_huber.wald gamma(1)=gamma(11)
freeze(wald_allbreak_ulc_pc)  model_orig_precovid_pcbreak_huber.wald gamma(2)=gamma(21)

freeze(wald_allbreak_nairu_var) model_orig_precovid_pcbreak_huber.wald @abs(sigma(10))=@abs(sigma(11))
freeze(wald_allbreak_nairu_var2) model_orig_precovid_pcbreak_huber.wald @abs(sigma(10)^2)=@abs(sigma(11)^2)
freeze(wald_allbreak_inf_var) model_orig_precovid_pcbreak_huber.wald @abs(sigma(1))=@abs(sigma(2))
freeze(wald_allbreak_inf_var2) model_orig_precovid_pcbreak_huber.wald @abs(sigma(1)^2)=@abs(sigma(2)^2)
freeze(wald_allbreak_ulc_var) model_orig_precovid_pcbreak_huber.wald @abs(sigma(3))=@abs(sigma(4))
freeze(wald_allbreak_ulc_var2) model_orig_precovid_pcbreak_huber.wald @abs(sigma(3)^2)=@abs(sigma(4)^2)

spool results_allbreaks

'Collating results for parameter breaks table coming from model_orig_precovid_pcbreak_huber
smpl 1968 2019
model_orig_precovid_pcbreak_huber.ml(m=500, cov=huber, showopts)

table(10,10)  param_breaks
param_breaks(1,2) = "Pre-break coef"
param_breaks(1,3) = "Post-break coef"
param_breaks(1,4) = "P value on break"


param_breaks(2,1)= "Inflation variance"
param_breaks(2,2)=sigma(1)^2
param_breaks(2,3)=sigma(2)^2
param_breaks(2,4)=wald_allbreak_inf_var2(6,4)

param_breaks(3,1)= "ULC variance"
param_breaks(3,2)=sigma(4)^2
param_breaks(3,3)=sigma(3)^2
param_breaks(3,4)=wald_allbreak_ulc_var2(6,4)

param_breaks(4,1)= "NAIRU variance"
param_breaks(4,2)=sigma(10)^2
param_breaks(4,3)=sigma(11)^2
param_breaks(4,4)=wald_allbreak_nairu_var2(6,4)

param_breaks(5,1)= "Inflation weight on pi expectations"
param_breaks(5,2)=(1-BETA(1)-BETA(2)-BETA(3)-BETA(6))
param_breaks(5,3)=(1-BETA(11)-BETA(21)-BETA(31)-BETA(61))
param_breaks(5,4)=wald_allbreak_inf_infx(6,4)

param_breaks(6,1)= "ULC weight on pi expectations"
param_breaks(6,2)=delta(4)
param_breaks(6,3)=delta(3)
param_breaks(6,4)=wald_allbreak_ulc_infx(6,4)

param_breaks(7,1)= "Inflation PC slope"
param_breaks(7,2)=gamma(1)
param_breaks(7,3)=gamma(11)
param_breaks(7,4)=wald_allbreak_inf_pc(6,4)


param_breaks(8,1)= "ULC PC slope"
param_breaks(8,2)=gamma(2)
param_breaks(8,3)=gamma(21)
param_breaks(8,4)=wald_allbreak_ulc_pc(6,4)

show param_breaks

 results_allbreaks.append model_orig_precovid_pcbreak_huber.results wald_allbreak_inf_infx wald_allbreak_ulc_infx wald_allbreak_inf_pc wald_allbreak_ulc_pc wald_allbreak_nairu_var wald_allbreak_inf_var wald_allbreak_ulc_var

show results_allbreaks

'Collate results for unrestrictred = 1 tests using model_orig_precovid_unres_huber
smpl 1968 2019
model_orig_precovid_unres_huber.ml(m=500, cov=huber, showopts)

' Run Wald test on unrestricted model (Technically F-tests better?)
freeze(wald_precovid_unres_infeq_postIT) model_orig_precovid_unres_huber.wald delta(1)+beta(1)+beta(2)=1  ' LR restriction on post-IT inflation eq
freeze(wald_precovid_unres_infeq_preIT) model_orig_precovid_unres_huber.wald delta(2)+beta(1)+beta(2)+beta(3)+beta(6)=1  ' LR restriction on pre-IT inflation eq
freeze(wald_precovid_unres_ulceq_postIT) model_orig_precovid_unres_huber.wald delta(3)=1  ' LR restriction on post-IT ULC eq
freeze(wald_precovid_unres_ulceq_preIT) model_orig_precovid_unres_huber.wald delta(4)+phi(4)=1  ' LR restriction on post-IT ULC eq

table(10,10) unrestricted_tests

unrestricted_tests(1,2)="Sum of lags + expectations"
unrestricted_tests(1,3)="P-value of restriction equal to 1"

unrestricted_tests(2,1)="Inflation pre-1993"
unrestricted_tests(2,2)=delta(2)+beta(1)+beta(2)+beta(3)+beta(6)
unrestricted_tests(2,3)=wald_precovid_unres_infeq_preIT(6,4)

unrestricted_tests(3,1)="Inflation post-1993"
unrestricted_tests(3,2)= delta(1)+beta(1)+beta(2) 
unrestricted_tests(3,3)=wald_precovid_unres_infeq_postIT(6,4)

unrestricted_tests(4,1)="ULC growth pre-1993"
unrestricted_tests(4,2)=  delta(4)+phi(4)
unrestricted_tests(4,3)=wald_precovid_unres_ulceq_preIT(6,4)

unrestricted_tests(5,1)="ULC growth post-1993"
unrestricted_tests(5,2)= delta(3)
unrestricted_tests(5,3)=wald_precovid_unres_ulceq_postIT(6,4)

show unrestricted_tests








spool results_breaks



param_breaks.save(t=xlsx, mode=update, cellfmt=eviews, strlen=256) "param_breaks.xlsx" RANGE=SHEET1!A1
unrestricted_tests.save(t=xlsx, mode=update, cellfmt=eviews, strlen=256) "unrestricted_tests.xlsx" RANGE=SHEET1!A1


