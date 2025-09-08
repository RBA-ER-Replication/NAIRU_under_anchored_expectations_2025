'===========================================================
' AER Policy Forum 2025: The NAIRU under anchored inflation expectations
' Main program for running NAIRU models and analysis for AER policy forum on NAIRU

' Ballantyne & Cusbert 2025 
'===========================================================
'-------------------------------------------------------------------------------------------------------
' Preliminaries and control
'-------------------------------------------------------------------------------------------------------

close @wf
mode quiet

%local_path = @runpath
cd %local_path

' Go up one level and into Data folder
%data_path = "..\Data\NAIRU_data.csv"

%out_path = %local_path + "\Output"
%prompt = "local (not shared)"
%data_start = "1964:1"

' Import NAIRU data, creates new workfile


import %data_path @freq q %data_start

string data_path = %data_path
string local_path = %local_path
string out_path = %out_path

cd %local_path


string data_start = %data_start
string prompt = %prompt



	' Dates for estimation etc, change as required
	string data_end    = u.@last         				' last period of data	
	string est_start  = "1968.1" 						' First period of estimation (1968:1 has become default)
	string est_end    = "2024.4" 						' Last period of estimation
	string IT_date = "1993:1" 								' Approximate introduction of the inflation target
	string NAV_date = "1984:1" 							' Nat Accounts GDP series break 
	string oil_date = "1976.4"								' Oil price dummy end date, seems arbitrary but as in Cusbert

	' Pandemic dates
	string covid_missing = "2020Q2 2020Q2"		' Dates for when you want missing data during height of pandemic
	string covid_ulc_extra = "2020Q4 2020Q4"		' Dates for when you want extra missing ULC data: 2021Q4 is borderline, residual is -3.5 which pulls NAIRU estimate down 0.2ppt
	string end_covid = "2023Q2"						' A date when you think the pandemic volatility is over
	string end_covid_short = "2022Q4"						' A date when you think the pandemic volatility is over
	string start_covid = "2020Q2" 

	
	' Which models do you want to run? Programs must follow existing naming convention: "model_[name].prg"
	string model_list = "orig_revised orig_precovid " ''orig_slope"  

	' Percentage trimming of start/end of estimation sample for rolling regs
	scalar trim_pct = 0.10    ' 0.05 to 0.50
	'string start_roll = "1980Q1"

' Various string
string wfname = "NAIRU_model_AERpaper" 	'Set name for workfile

' Cut workfile at end of sample in case trailing data
pagestruct(end={est_end})

' Save workfile, active directory is @runpath
%destwf = wfname+".wf1"
wfsave %destwf 



'-------------------------------------------------------------------------------------------------------
' Data transformations
'-------------------------------------------------------------------------------------------------------


' Standard transformations
smpl @all
' Quarterly

series pmq=pm/pm(-1)*100-100
series infq=inf/inf(-1)*100-100
series infxq=infx/4
series ulcq=ulc/ulc(-1)*100-100
' Year ended
series infye=inf/inf(-4)*100-100
series pmye=pm/pm(-4)*100-100
series ulcye=ulc/ulc(-4)*100-100

' Inflation targeting dummy
smpl @all
series IT=0
smpl {IT_date} @last
IT=1

' National Accounts Volatilty: New GDP series break
smpl @all
series NAV=1
smpl {NAV_date} @last
NAV=0

' Oil dummy prior to 1977
smpl @all
series oil_dummy=0
smpl @first {oil_date}
oil_dummy=1


' Make core pandemic missing data, suprisingly only need a few quarters for outliers
smpl @all
copy infq infq_orig
'copy aenaq aenaq_orig
copy ulcq ulcq_orig
'copy wpiq wpiq_orig

smpl {covid_missing}
'infq = NA
'aenaq = NA
ulcq = NA
'wpiq = NA
' Additional ULC drop
smpl {covid_ulc_extra}
ulcq = NA
smpl @all


' Make volatility break -
smpl @all
series COVID = 0
smpl 2020Q2 {end_covid}
series COVID = 1
smpl @all

' Make infq error persistence break - from when inflation starts to pick up due to cost-push
smpl @all
series COVID_inf_error = 0
smpl 2021Q2 {end_covid} 
series COVID_inf_error  = 1
smpl @all



'--------------------------------------------------------------------------------
' Call model routines
'--------------------------------------------------------------------------------

' Loop over models 
for %m {model_list}

	' Run model
	%temp = local_path+"\model_"+%m+".prg"
	exec %temp

	' Do other stuff - future infrastructure
	
next

%temp1 = local_path+"Uncertainty/nairu_uncertainty_main.prg"
exec %temp1



'--------------------------------------------------------------------------------
' Test biased expectations
'--------------------------------------------------------------------------------

' Copy model, stash pie and shift
copy model_orig model_orig_pieshift
smpl {est_start} {est_end}
series infxq_orig = infxq 
series infxq = infxq +(0.5/4)

' Run bias pie model
include check_sspace.prg ' Include sspace check
model_orig_pieshift.ml(m=500, showopts)
call check_sspace("model_orig_pieshift")
' Create smoothed nairu and standard errors 
model_orig_pieshift.makestates(t=smooth) *_orig_pieshift_sm
'line nairu_orig_sm nairu_orig_pieshift_sm

' Reset pie, output
series infxq = infxq_orig
series nairu_diff_pieshift = nairu_orig_pieshift_sm - nairu_orig_sm
nairu_diff_pieshift.stats 

'--------------------------------------------------------------------------------
' Run scenarios anchored/adaptive inflation & make PC scatter data
'--------------------------------------------------------------------------------

' Run anchored_decay  after running this script because it overwrites u data

'--------------------------------------------------------------------------------
' Extras
'--------------------------------------------------------------------------------

'------------------------------------------------
' Graph data for COVID residuals graph
' Fixed parameters on precov model, restore original ULCs
smpl @all
model_orig_precovid.makefilter filter_orig_precovid
series ulcq_est = ulcq
ulcq = ulcq_orig
' Run fixed param model over full sample, get resids
smpl {est_start} {est_end}
filter_orig_precovid.ml(m=500,showopts)
filter_orig_precovid.makesignals(t=resid) *_resid_precov
' Revert data
smpl @all
ulcq = ulcq_est

smpl 2000 2025
show date infq_resid_precov date ulcq_resid_precov

'code for graphs about volatility
smpl 1968 2019
GROUP graphdata_g4g5 date infye infx date infq_resid_precov  date ulcye infx  date ulcq_resid_precov

%temp_output = local_path +"Output\inf_ulc_graphdatag4_5.xlsx"
' Write scatter data to excel

wfsave(type=excelxml, noid) %temp_output @keep dates  graphdata_g4g5 @smpl 1968 2019

graph g4a_inflationvol infye infx
show g4a_inflationvol

graph g4b_inflationresid infq_resid_precov  
show g4b_inflationresid 

graph g5a_ulcvol ulcye infx
show g5a_ulcvol

graph g5b_ulcresid ulcq_resid_precov
show g5b_ulcresid 

smpl 2000 2024
graph g11a_inflationcovid  infq_resid_precov
show g11a_inflationcovid

graph g11b_ulccovid  ulcq_resid_precov
show g11b_ulccovid

'------------------------------------------------
' Graph data for Estimates of Supply-driven inflation
' Make AR component of inf error terms
smpl @all
series ar_e_est = COVID_INF_ERROR*coef_orig_zeta(99)*ar_e_orig_sm(-1)
series ar_e_est_ye = (1+ar_e_est/100)*(1+ar_e_est(-1)/100)*(1+ar_e_est(-2)/100)*(1+ar_e_est(-3)/100)*100-100
series ar_e_orig_smye = (1+ar_e_orig_sm/100)*(1+ar_e_orig_sm(-1)/100)*(1+ar_e_orig_sm(-2)/100)*(1+ar_e_orig_sm(-3)/100)*100-100

smpl 2020 2025
'show infye ar_e_orig_ye 
show infye ar_e_est_ye ar_e_orig_smye
'show infq ar_e_est ar_e_orig_sm


'-------------------------------------------------------------------------------------------------------
' Save and exit stage left
'-------------------------------------------------------------------------------------------------------


'Appedinx table

freeze(Appendixtable) model_orig_huber
show appendixtable

%save_location = local_path
cd %save_location


' Always save workfile in active directory @runpath
%destwf = wfname+".wf1"
wfsave %destwf

'Run anchored decay scenarios last

  %temp2 = local_path+"anchored_decay.prg"
exec %temp2

' NB below version has full set of results but incorrect u data
%destwf = wfname+"_scenarios.wf1"
wfsave %destwf


