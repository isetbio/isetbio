function [fit_out,x,err] = FitGammaExtP(values_in,measurements,values_out,x0)
% [fit_out,x,err] = FitGammaExtP(values_in,measurements,values_out,x0)
%
% Fit extended power function to gamma data.
%
% 4/08/02 awi   Turned off warnings before calling constr to avoid future obsolete warning.  
% 4/13/02 dgp   Fixed the warning suppression to save and restore original state.
% 3/4/05  dhb   Handle new version of optimization toolbox.
% 8/2/15  dhb   Update calls to exist to work when query is for a p-file.

% Check for needed optimization toolbox, and version.
if (exist('fminunc','file'))
    options = optimset;
    options = optimset(options,'Display','off');
    if ~IsOctave
        options = optimset(options,'LargeScale','off');
    end

    x = fminunc(@FitGammaExtPFun,x0,options);
elseif (exist('constr','file'))
    options = foptions;
    options(1) = 0;
    options(14) = 600;
    x = constr(@FitGammaExtPFun,x0,options);
else
    error('FitGammaExtP requires the Matlab Optimization Toolbox from Mathworks');
end

% Now compute fit values and error to data for return
fit_out = ComputeGammaExtP(x,values_out);
err = FitGammaExtPFun(x);

    % Nested target function to optimize:
    function [err,con] = FitGammaExtPFun(x)
        % [err,con] = FitGammaExtPFun(x)
        % 
        % Error function for power function fit.
        %
        % 9/21/93  dhb  Added positivity constraint.
        % 8/12/15  dhb  Fix typo, ComputeFSWSE -> ComputeFFSE.

        predict = ComputeGammaExtP(x,values_in);
        err = ComputeFSSE(measurements,predict);
        con = [-x];
    end
end
