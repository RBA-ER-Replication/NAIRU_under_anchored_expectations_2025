'===========================================================
' Estimating filtering and parameter uncertainty for NAIRU models
' Main program - built so that you can rerun the model blocks easily
' Ballantyne Jan 2024, using Matt Read's uncertainty code
'===========================================================

'-------------------------------------------------------------------------------------------------------
' Preliminaries and control
'-------------------------------------------------------------------------------------------------------



%path = @runpath
cd %path


cd "Uncertainty"

%path = %path+"Uncertainty"



' Switch workfile if needed
'if %shared_wf = "1" then
'	%wfswitch = shared_path+"\Output\"+wfname
'	wfclose 
'	wfopen %wfswitch
'endif

'-------------
' IMPORTANT
' Specify the measure of slack, 
%v = "u" 'Measure of labour market slack in Phillips curves: 'u = unemployment rate or util = underutilisation rate
' Do you want a partial estimate one quarter after the most recent set of full data?
string partial = "no"   ' yes or no
'-------------

' Control strings
string uncertdir = %path
string graphdir = out_path
string outputuncert = "Uncertainty_output.xlsx" 'Set name for Excel output
string wfnameuncert = "NAIRU_Uncertainty_"+%v 	'Set name for workfile when saving
string slack_measure = %v  ' to pass to modules
scalar shared_wf = @val(%shared_wf)

' Which models do you want to run? Programs must follow existing naming convention: "model_[name].prg"
string uncert_list = "orig"

'Specify number of Monte Carlo replications used to approximate MSE of state estimates
scalar K = 1000
rndseed 23032021

' Make better date to export
alpha dates = @datestr(@dateadd(@date,2,"mm"), "DD/MM/YYYY")

'Save workfile to change name
%destwf = uncertdir+"\"+wfnameuncert+".wf1"
wfsave %destwf

' For later use
%basewf = wfnameuncert

' Cut workfile at end of sample 
pagestruct(end={est_end})


'--------------------------------------------------------------------------------
' Call model routines
'--------------------------------------------------------------------------------

' Loop over NAIRU models 
for %m {uncert_list}

	'Specify model variant (e.g. wpi) and adjust for slack measure
	string varname = %m
	if %v = "util" then
		varname = varname+"_"+%v
	endif
	%varname = varname	

	' Run model
	%temp = uncertdir+"\uncertainty_"+%m+"_revised.prg"
	exec %temp

	' Some distribution output, mostly for note
	scalar rate = @elem({slack_measure},{slack_measure}.@last)
	wfcreate(wf=unstruct, page=undated) u {k} ' create workfile with K obs
	copy {%basewf}::Nairu_data\probs_{%varname} unstruct::undated\probs_{%varname}  ' Copy required objects
	copy {%basewf}::Nairu_data\draws_{%varname}_sm unstruct::undated\draws_{%varname}_sm
	copy {%basewf}::Nairu_data\rate unstruct::undated\rate  
	vector currentp_{%varname} = probs_{%varname}.@row(probs_{%varname}.@rows) ' select last row
	vector currentn_{%varname} = draws_{%varname}_sm.@row(draws_{%varname}_sm.@rows) ' select last row
	mtos(currentp_{%varname}, currentp_{%varname}_series) ' turn into series
	mtos(currentn_{%varname}, currentn_{%varname}_series)
	series currentg_{%varname}_series = rate - currentn_{%varname}_series ' make gap instead of nairu
	currentp_{%varname}_series.distdata(dtype=hist, anchor=0, binw=user, binval=0.01) histp_{%varname}  ' make histogram data
	currentg_{%varname}_series.distdata(dtype=hist, anchor=0, binw=user, binval=0.1) histg_{%varname}
	copy unstruct::undated\histp_{%varname} {%basewf}::Nairu_data\hist_prob_{%varname} 
	copy unstruct::undated\histg_{%varname} {%basewf}::Nairu_data\hist_gap_{%varname} 
	wfclose unstruct
 
next


'----------------------------------
' SAVE OUTPUT - sometimes it doesn't write to excel properly so can just run this segment

' Remake dates over full range
smpl @all
alpha dates = @datestr(@dateadd(@date,2,"mm"), "DD/MM/YYYY")


	' Loop to save models output
	for %m {uncert_list}
	
		'Specify model variant (e.g. wpi) and adjust for slack measure
		string varname = %m
		%sheetname = varname+"_"+%v  'Won't write three letter sheet names would you believe! 
		if %v = "util" then
			varname = varname+"_"+%v
		endif
		%varname = varname	

smpl 1980 2024 	
series SE_total_orig_smooth= (smoother_var_orig+par_var_orig)^0.5
group graph10data u nairu_orig_sm  nairu_orig_sm+1*SE_total_orig_smooth nairu_orig_sm+2*SE_total_orig_smooth nairu_orig_sm-1*SE_total_orig_smooth nairu_orig_sm-2*SE_total_orig_smooth

graph g10a_uncertainty_nairu graph10data
show g10a_uncertainty_nairu

graph g10b_prob_nairu probs_avg_orig
show g10b_prob_nairu
		'Export results to Excel

	%destout = graphdir+"\uncertainty_graph_data"
wfsave(type=excelxml, mode=update, noid) %destout range = sheet1 @keep date graph10data probs_avg_orig @smpl  1980 2024

		  '  +%varname+%v+".xlsx"
			group results {%v} nairu_{%varname}_sm nairu_{%varname}_smse smoother_var_{%varname} par_var_{%varname} probs_avg_{%varname}
			'write(t=xls, name) %destout results
			if partial = "yes" then
				
			else
				
			endif
	
	next



'Save workfile locally always
'%destwf = uncertdir+"\"+wfnameuncert+".wf1"
'wfsave %destwf


