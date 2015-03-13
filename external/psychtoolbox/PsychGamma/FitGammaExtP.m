function [fit_out,x,err] = FitGammaPow(values_in,measurements,values_out,x0)
% [fit_out,x,err] = FitGammaPow(values_in,measurements,values_out,x0)
%
% Fit extended power function to gamma data.
%
% 4/08/02 awi	Turned off warnings before calling constr to avoid future obsolete warning.  
% 4/13/02 dgp	Fixed the warning suppression to save and restore original state.
% 3/4/05  dhb   Handle new version of optimization toolbox.

% Check for needed optimization toolbox, and version.
if (exist('fmincon') == 2)
	options = optimset;
	options = optimset(options,'Diagnostics','off','Display','off');
	options = optimset(options,'LargeScale','off');
	x = fminunc('FitGammaExtPFun',x0,options,values_in,measurements);	
elseif (exist('constr') == 2)
	options = foptions;
	options(1) = 0;
	options(14) = 600;
	x = constr('FitGammaExtPFun',x0,options,[],[],[],values_in,measurements);
else
	error('FitGammaExtP requires the Matlab Optimization Toolbox from Mathworks');
end

% Now compute fit values and error to data for return
fit_out = ComputeGammaExtP(x,values_out);
err = FitGammaExtPFun(x,values_in,measurements);
