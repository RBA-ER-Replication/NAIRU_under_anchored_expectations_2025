'===========================================================
' Work for AER Policy Forum on NAIRU 2025
' Show what pi looks like under different unemployment rate scenarios in short-run and under adaptive expectations.

' Ballantyne/Cusbert Jan 2025
'===========================================================

subroutine make_ye (series s)

' Make ye series
%name = s.@name
series {%name}_ye = (1+s/100)*(1+s(-1)/100)*(1+s(-2)/100)*(1+s(-3)/100)*100-100

endsub


'--------------------------------------------------------------------------------
' 1. Setup 
'--------------------------------------------------------------------------------

' Extend range
pagestruct(end=2030Q4)

' Settings ----------------------------------------------------

' Model name
%v = "u"
%mname = "orig"

' Start of scenario (actually last period of actual data)
string scn_start = "2022:4"

' Peak unemployment in overshooting scenario
scalar u_peak = 5.75

'------------------------------------------------------------------

' Set up parameters for ease
c = coef_{%mname}_c
beta = coef_{%mname}_beta
delta = coef_{%mname}_delta
gamma = coef_{%mname}_gamma
omega = coef_{%mname}_omega
psi = coef_{%mname}_psi
phi = coef_{%mname}_phi

' State initialization vectors
vector sprior_scn = sprior_{%mname}
sym vprior_scn = vprior_{%mname}
sprior_scn(1) = @elem(nairu_{%mname}_sm,scn_start)
sprior_scn(2) = 0
vprior_scn(1,1) = @elem(nairu_{%mname}_smse,scn_start)^2
vprior_scn(2,2) = 0


' Last data
string samp_end = infq.@last

' Some maths
scalar pibar = 2.5
scalar pibarq = (1+pibar/100)^(1/4)*100-100
scalar pibarqdl = log(1+pibarq/100)

' Set rates | exog variables set to convenient theoretical steady state
scalar nairu_baseline = @elem(nairu_{%mname}_sm,scn_start)
scalar nairu_high = nairu_baseline+0.5
scalar nairu_low = nairu_baseline-0.5
scalar nairu_central = 4.5
scalar u_baseline = @elem(u,scn_start)
scalar u_high = u_baseline+1
scalar u_low = u_baseline-1
scalar infq_baseline = @elem(infq,scn_start)
scalar infq_high = infq_baseline+1
scalar infq_low = infq_baseline-1
scalar infx_baseline = pibar ' @elem(infx,infx.@last)
scalar infxq_baseline = pibarq ' @elem(infxq,infxq.@last)
scalar pmye_baseline = infx_baseline ' @elem(pmye,pmye.@last)
scalar oil_baseline = @elem(oil,scn_start)


' Extend dummies
smpl {samp_end} @last
COVID = 0
NAV = 0
IT = 1
OIL_DUMMY = 0
covid_inf_error = 0

' Remake dates over full range
smpl @all
alpha dates = @datestr(@dateadd(@date,2,"mm"), "DD/MM/YYYY")
series target = pibar


' Save orig data then extend
smpl @all
series u_old = u
series infq_old = infq
series infx_old = infx
series infxq_old = infxq
'series laq_old = laq
'series aenaq_old = aenaq
series ulcq_old = ulcq
series pmye_old = pmye
series oil_old = oil
'series wpiq_old = wpiq
'series py_old = py
series nairu = nairu_orig_sm  'make orig
series nairu_old = nairu

' Make baseline paths for RHS variables other than u
smpl {scn_start}+1 @last
series infx = infx_baseline
series infxq = infxq_baseline
series pmye = pmye_baseline
series oil = oil_baseline


'--------------------------------------------------------------------------------
' 2. Project - anchored expectations
'--------------------------------------------------------------------------------

' Inflation equation for reference
'@SIGNAL INFQ = IT*(1-BETA(1)-BETA(2))*INFXQ+(1-IT)*(1-BETA(1)-BETA(2)-BETA(3)-BETA(6))*INFXQ + BETA(1)*INFQ(-1) + BETA(2)*INFQ(-2) + (1-IT)*BETA(3)*INFQ(-3) + (1-IT)*BETA(6)*(AENAQ(-1)-TREND_PROD_LAG) + GAMMA(1)*(U-NAIRU)/U + (1-IT)*OMEGA(1)*D(U(-1))/U + PSI(1)*(PMYE(-1) - INFX(-1)) + C(2)*OIL_DUMMY*D(LOG(OIL(-2))) + AR_E

'---------------------------------------
' Baseline 

' Make unemployment path
smpl @first {scn_start}
series u_const = u_old
smpl {scn_start}+1 @last ' Hold u constant
u_const = u_baseline
' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_const

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' Ridiculous Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_baseline_* @signal proj_baseline_* 
	for %s infq ulcq nairu ' would also need to do states if other than random walk
		{%s} = proj_baseline_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq nairu
	proj_baseline_{%s} = {%s}
	{%s} = {%s}_old
next



'---------------------------------------
' Unemp @ NAIRU

' Make unemployment path
smpl @first {scn_start}
series u_nairu = u_old
smpl {scn_start}+8 @last ' Return to NAIRU 8 quarters after scenario
u_nairu = nairu_baseline
'Interpolate gaps
smpl {scn_start} @last 
u_nairu.ipolate(type=cb) u_nairu_int
u_nairu = u_nairu_int
' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_nairu

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' NB Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_nairu_* @signal proj_nairu_* 
	for %s infq ulcq nairu ' would also need to do states if other than random walk
		{%s} = proj_nairu_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq  nairu
	proj_nairu_{%s} = {%s}
	{%s} = {%s}_old
next


'---------------------------------------
' Unemp overshoot scenario

' Make unemployment path
smpl @first {scn_start}
series u_scn = u_old
smpl {scn_start}+4 {scn_start}+4 ' Overshoot 4 quarters after scenario
u_scn = u_peak
smpl {scn_start}+8 @last ' Return to NAIRU 8 quarters after scenario
u_scn = nairu_baseline
'Interpolate gaps
smpl {scn_start} @last 
u_scn.ipolate(type=cb) u_scn_int
u_scn = u_scn_int
' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_scn

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' NB Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_scn_* @signal proj_scn_* 
	for %s infq ulcq nairu ' would also need to do states if other than random walk
		{%s} = proj_scn_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq  nairu
	proj_scn_{%s} = {%s}
	{%s} = {%s}_old
next


'--------------------------------------------------------------------------------
' 3. Project - adaptive expectations
'--------------------------------------------------------------------------------


'---------------------------------------
' Baseline with adaptive expectations

' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_const

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' Ridiculous Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	if !i > 4 then
		series infxq = @movav(infq(-1),4)
	endif
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_adaptbase_* @signal proj_adaptbase_* 
	for %s infq ulcq  nairu ' would also need to do states if other than random walk
		{%s} = proj_adaptbase_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq nairu
	proj_adaptbase_{%s} = {%s}
	{%s} = {%s}_old
next

'---------------------------------------
' Unemp @ NAIRU with adaptive expectations

' Make unemployment path
' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_nairu

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' NB Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	if !i > 4 then
		series infxq = @movav(infq(-1),4)
	endif
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_adaptnairu_* @signal proj_adaptnairu_* 
	for %s infq ulcq nairu ' would also need to do states if other than random walk
		{%s} = proj_adaptnairu_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq  nairu
	proj_adaptnairu_{%s} = {%s}
	{%s} = {%s}_old
next

'---------------------------------------
' Unemp overshoot scenario with adaptive expectations

' Overwrite u variable to feed model
smpl {scn_start}+1 @last
series u = u_scn

smpl {scn_start}+1 @last
for !i = 1 to @obssmpl  ' Ridiculous Eviews cannot do dynamic state space forecasting properly
	smpl {scn_start}+!i {scn_start}+!i
	if !i > 4 then
		series infxq = @movav(infq(-1),4)
	endif
	' Forecasting can use i=o option if after the estimation sample, or i=u option if before
	model_{%mname}.forecast(i=u, m=d, mprior=sprior_scn, vprior=vprior_scn) @state proj_adaptscn_* @signal proj_adaptscn_* 
	for %s infq ulcq nairu ' would also need to do states if other than random walk
		{%s} = proj_adaptscn_{%s}
	next
next
' Stash and reset
smpl @all
for %s infq ulcq  nairu
	proj_adaptscn_{%s} = {%s}
	{%s} = {%s}_old
next


'--------------------------------------------------------------------------------
' 4. Viz
'--------------------------------------------------------------------------------

'Splice series
smpl @all
for %s inf ulc 
	series {%s}_baseline = {%s}q
	series {%s}_nairu = {%s}q
	series {%s}_scn = {%s}q
	series {%s}_adaptbase = {%s}q
	series {%s}_adaptnairu = {%s}q
	series {%s}_adaptscn = {%s}q
next
smpl {scn_start}+1 @last
for %s inf ulc 
	{%s}_baseline = proj_baseline_{%s}q
	{%s}_nairu = proj_nairu_{%s}q
	{%s}_scn = proj_scn_{%s}q
	{%s}_adaptbase = proj_adaptbase_{%s}q
	{%s}_adaptnairu = proj_adaptnairu_{%s}q
	{%s}_adaptscn = proj_adaptscn_{%s}q
next

' Make year-ended
smpl @all
for %s inf_baseline inf_nairu inf_scn ulc_baseline ulc_nairu ulc_scn 
	call make_ye ({%s})
next
smpl @all
for %s inf_adaptbase inf_adaptnairu inf_adaptscn ulc_adaptbase ulc_adaptnairu ulc_adaptscn 
	call make_ye ({%s})
next



smpl 2020 2027

graph g7_uratescens.line  u_const u_nairu u_scn 
show g7_uratescens

graph g8_infanch.line  inf_baseline_ye inf_nairu_ye inf_scn_ye 
show  g8_infanch

graph g9_infadapt.line  inf_adaptbase_ye inf_adaptnairu_ye inf_adaptscn_ye 
show g9_infadapt

'--------------------------------------------------------------------------------
' 5. Make Phillips curve
'--------------------------------------------------------------------------------

' Make u gap
smpl @all
series gap_orig = u-nairu_orig_sm

' Create simulated values for unemployment and pi 
smpl 1993Q1 1993Q1+65
delete(noerr) ugap_sim_vec 
vector u_sim_vec = @seq(3, 0.1, @obssmpl) 	
mtos(u_sim_vec, u_sim)

' Fixed 
series ugap_sim = u_sim-nairu_central
series pi_sim = 2.5+COEF_orig_GAMMA(1)*4*(ugap_sim)/u_sim ' annualise with times 4, makes minimal difference doing it exact compounding
scalar pie_coef_orig =1-COEF_orig_BETA(1)-COEF_orig_BETA(2)
series pi_sim_lr = 2.5+COEF_orig_GAMMA(1)*(pie_coef_orig)*4*(ugap_sim)/u_sim ' Cusbert long-run version 
smpl 1993Q1 {est_end}


' Rinse and repeat for full history
smpl {est_start} {est_end}
delete(noerr) ugap_sim_vec_full
vector u_sim_vec_full = @seq(1.5, 9/@obssmpl, @obssmpl) 	
mtos(u_sim_vec_full, u_sim_full)
' Fixed 
series ugap_sim_full = u_sim_full-nairu_central
series pi_sim_full = 2.5+COEF_orig_GAMMA(1)*4*(ugap_sim_full)/u_sim_full ' annualise with times 4, makes minimal difference doing it exact compounding
series pigap_sim_full = COEF_orig_GAMMA(1)*4*(ugap_sim_full)/u_sim_full
series infgap = infye-infx
graph g1_infPC.scatpair ugap_sim_full pigap_sim_full gap_orig infgap ' ugap_sim pi_sim_lr
show g1_infPC
' Rince and repeat for ULC
smpl {est_start} {est_end}
call make_ye (ulcq_orig)
series ulcgap_sim_full = COEF_ORIG_GAMMA(2)*4*(ugap_sim_full)/u_sim_full
series ulcgap = ulcq_orig_ye-infx
graph g2_ULCPC.scatpair ugap_sim_full ulcgap_sim_full gap_orig ulcgap ' ugap_sim pi_sim_lr
show  g2_ULCPC

' Show data fro Graphit
smpl {est_start} {est_end}
show dates gap_orig infgap ugap_sim_full pigap_sim_full gap_orig ulcgap ugap_sim_full ulcgap_sim_full  




'--------------------------------------------------------------------------------
' 5. Create outputs
'--------------------------------------------------------------------------------

'smpl 2020 2027
'show dates inf_baseline_ye inf_nairu_ye inf_adaptbase_ye infye

'show inf_nairu_ye aena_nairu_ye wpi_nairu_ye

'smpl 1993Q1 {est_end}
'show gap_orig infye ugap_sim pi_sim


%temp_output = local_path +"Output\anchored_decay.xlsx"
' Write anchor/adaptive scenarios data to excel
group scn_results u_old u_const u_nairu u_scn infye inf_baseline_ye inf_nairu_ye inf_scn_ye inf_adaptbase_ye inf_adaptnairu_ye inf_adaptscn_ye ulcye ulc_baseline_ye ulc_nairu_ye ulc_scn_ye ulc_adaptbase_ye ulc_adaptnairu_ye ulc_adaptscn_ye 
wfsave(type=excelxml, noid) %temp_output  @keep dates scn_results @smpl 2015 2027

%temp_output = local_path +"Output\PCscatters.xlsx"
' Write scatter data to excel
group scat_results gap_orig infgap ugap_sim_full pigap_sim_full gap_orig ulcgap ugap_sim_full ulcgap_sim_full 
wfsave(type=excelxml, noid) %temp_output @keep dates scat_results @smpl {est_start} {est_end}


