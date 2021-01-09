function [fit_out,x,err] = FitGammaSig(values_in,measurements,values_out,x0)
% [fit_out,x,err] = FitGammaSig(values_in,measurements,values_out,x0)
%
% Fit sigmoid function to gamma data.
%
% 3/4/05    dhb Conditionals for optimization toolbox version.
% 8/2/15    dhb Update calls to exist to work when query is for a p-file.

% Check for needed optimization toolbox, and version.
if (exist('fminunc','file'))
    options = optimset;
    options = optimset(options,'Display','off');
    if ~IsOctave
        options = optimset(options,'LargeScale','off');
    end

    x = fminunc(@FitGammaSigFun,x0,options);
elseif (exist('constr','file'))
    options = foptions;
    options(1) = 0;
    options(14) = 600;
    x = constr(@FitGammaSigFun,x0,options);
else
    error('FitGammaSig requires the optional Matlab Optimization Toolbox from Mathworks');
end

% Now compute fit values and error to data for return
fit_out = ComputeGammaSig(x,values_out);
err = FitGammaSigFun(x);

    % Nested target function to optimize:
    function [err,con] = FitGammaSigFun(x)
    % [err,con] = FitGammaSigFun(x)
    % 
    % Error function for sigmoid function fit.

    predict = ComputeGammaSig(x,values_in);
    err = ComputeFSSE(measurements,predict);
    con = [-x(1) -x(2) -x(3)];
    end
end
