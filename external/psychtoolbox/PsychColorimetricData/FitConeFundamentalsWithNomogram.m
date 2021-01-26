function [params,fitFundamentals,fitError] = FitConeFundamentalsWithNomogram(T_targetQuantal,staticParams,params0)
% [fitFundamentals,params,fitError] = FitConeFundamentalsWithNomogram(T_targetQuantal,staticParams,params0)
%
% Find underlying parameters that fit the passed corneal cone fundamentals.
%
% Needs the Matlab optimization toolbox, or GNU/Octave version 6 or later with
% the Octave Forge 'optim' package v1.6.0 or later.
%

% 8/4/03   dhb  Wrote it.
% 12/4/20  mk  Add Octave support.

% On Octave we need the 'optim' package v1.6.0 or later for fmincon() support,
% and Octave 6 or later for support for handles to nested functions like @FitConesFun below:
if IsOctave
    v = version;
    if str2num(v(1)) < 6
        error('For use with Octave, you need at least Octave version 6.');
    end

    try
        % Try loading the optim package with the optimization functions:
        pkg load optim;
    catch
        error('For use with Octave, you must install the ''optim'' package from Octave Forge. See ''help pkg''.');
    end

    % Got optim package loaded. Does it support fmincon()?
    if ~exist('fmincon')
        error('For use with Octave, you need at least version 1.6.0 of the ''optim'' package from Octave Forge.');
    end
else
    if ~exist('fmincon')
        error('For use with Matlab, you need the optimization toolbox.');
    end
end

% Convert initial parameter struct to parameter list
x0 = FitConesParamsToList(params0);

% Set bounds on search parameters
% Length 4, we're handling Ser/Ala polymorphism
if (length(x0) == 4)
    vlb = [0 ; 0; 0; 0];
    vub = [800; 800; 800; 800];
elseif (length(x0) == 3)
    vlb = [0 ; 0; 0];
    vub = [800; 800; 800];
else
    error('Unexpected length for parameter vector');
end

% Search to find best fit
options = optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','active-set');
if (exist('IsCluster') && IsCluster && matlabpool('size') > 1) %#ok<EXIST>
    options = optimset(options,'UseParallel','always');
end
x = fmincon(@FitConesFun,x0,[],[],[],[],vlb,vub,[],options);

% Convert parameter list to parameter struct and
% compute final values for return
params = FitConesListToParams(x);
fitError = FitConesFun(x);
fitFundamentals = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
            staticParams.pupilDiameterMM,params.lambdaMax,staticParams.whichNomogram ...
            );

    function [f] = FitConesFun(x)
        DO_LOG = 1;
        params = FitConesListToParams(x);
        T_pred = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
            staticParams.pupilDiameterMM,params.lambdaMax,staticParams.whichNomogram ...
            );
        
        if (DO_LOG)
            bigWeight = 1; bigThresh = -1;
            index = find(~isinf(log10(T_targetQuantal(:))));
            T_resid = log10(T_pred(index))-log10(T_targetQuantal(index));
            index1 = log10(T_targetQuantal(index)) > bigThresh;
            index2 = log10(T_targetQuantal(index)) <= bigThresh;
            if ( any(isnan(T_resid(:))) || any(isinf(T_resid(:))) )
                f = 1e6;
            else
                f = 100*(bigWeight*mean(T_resid(index1).^2) + mean(T_resid(index2).^2));
            end
        else
            T_resid = T_pred-T_targetQuantal;
            f = 100*mean((T_resid(:)).^2);
        end
    end
end

% Convert list to parameter structure
function params = FitConesListToParams(x)

% Length 4, we're handling Ser/Ala polymorphism
if (length(x) == 4)
    params.lambdaMax = x(1:4);
elseif (length(x) == 3)
    params.lambdaMax = x(1:3);
else
    error('Unexpected length for parameter vector');
end

end

% Convert parameter structure to list
function x = FitConesParamsToList(params)

if (size(params.lambdaMax,1) == 4)
    x = zeros(4,1);
    x(1:4) = params.lambdaMax;
elseif (size(params.lambdaMax,1) == 3)
    x = zeros(3,1);
    x(1:3) = params.lambdaMax;
else
    error('Unexpected number of photopigment lambda max values passed');
end

end
