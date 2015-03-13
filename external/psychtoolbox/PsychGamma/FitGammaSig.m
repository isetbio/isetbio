function [fit_out,x,err] = FitGammaSig(values_in,measurements,values_out,x0)
% [fit_out,x,err] = FitGammaSig(values_in,measurements,values_out,x0)
%
% Fit sigmoid function to gamma data.

% 3/4/05	dhb	  Conditionals for optimization toolbox version.

% Check for needed optimization toolbox, and version.
if (exist('fmincon') == 2)
	options = optimset;
	options = optimset(options,'Diagnostics','off','Display','off');
	options = optimset(options,'LargeScale','off');
	x = fminunc('FitGammaSigFun',x0,options,values_in,measurements);	
elseif (exist('constr') == 2)
	options = foptions;
	options(1) = 0;
	options(14) = 600;
	x = constr('FitGammaSigFun',x0,options,[],[],[],values_in,measurements);
else
	error('FitGammaSig requires the optional Matlab Optimization Toolbox from Mathworks');
end

% Now compute fit values and error to data for return
fit_out = ComputeGammaSig(x,values_out);
err = FitGammaSigFun(x,values_in,measurements);
