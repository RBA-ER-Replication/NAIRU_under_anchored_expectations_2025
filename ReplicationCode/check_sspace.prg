'--------------------------------------------------------------------------------
' Subroutine for sspace convergence checks, Ballantyne Feb 2024
'--------------------------------------------------------------------------------
subroutine check_sspace (string %modname)

'Check estimation convergence, a few ways
vector zstats = {%modname}.@tstats  ' get z-stats
vector zbad = @egt(@abs(zstats),100*@ones(zstats.@rows))  ' z-stats > 100 get flagged
vector zna = @eeqna(zstats,na*@ones(zstats.@rows))  ' z-stats that are nas trigger exit
scalar zcount = @sum(zbad)
scalar zcountna = @sum(zna)
vector pvals = 2*@cnorm(-@abs(zstats))  ' calc p-values
vector pbad = @egt(@abs(pvals),0.9999*@ones(zstats.@rows))  ' p-vals > 0.999 get flagged
scalar pcount = @sum(pbad)
freeze(sscheck) {%modname}.output
string sconv = sscheck.@find("@instr([@all],""Convergence not achieved"")") ' check for running out of iterations
delete sscheck 
%ssname = {%modname}.@displayname  ' name of model
if sconv <> "" then
	@uiprompt(%ssname+" state space model did not converge to a solution. Please check the model and diagnose the issue. Exiting program.", "O")
	stop
else 
if zcountna > 0 then
	@uiprompt(%ssname+" state space model has z-stats with NA values so the model did not converge. Please check the model and diagnose the issue. Exiting program.", "O")
	stop
else
if pcount > 0 then
	@uiprompt(%ssname+" state space model has p-values greater than 0.999, which is a serious concern. Please check the model and diagnose the issue. Exiting program.", "O")
	stop
else
if zcount > 0 then
	@uiprompt(%ssname+" state space model has z-stats greater than 100, which is a serious concern. Please check the model and diagnose the issue. Exiting program.", "O")
	stop
endif
endif
endif
endif

endsub


