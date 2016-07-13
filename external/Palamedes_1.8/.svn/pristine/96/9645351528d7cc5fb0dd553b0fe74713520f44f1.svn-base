%
%PAL_PFLR_TLR       Fit two models to a multi-condition experiment and
%   determine the value of the transformed likelihood ratio.
%
%syntax: [TLR exitflag paramsL paramsF] = PAL_PFLR_TLR(StimLevels, ...
%   NumPos, OutOfNum, params, PF, MC, {optional arguments})
%
%Internal function
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.1.0, 1 .3.0, 1.4.0, 1.6.3 (see History.m)

function [TLR, exitflag, paramsL, paramsF, funcParamsL, funcParamsF] = PAL_PFLR_TLR(StimLevels, NumPos, OutOfNum, paramsValues, PF, MC, varargin)

if size(paramsValues,1) == 1
    paramsValues = repmat(paramsValues,[size(StimLevels,1) 1]);
end

alphas = paramsValues(:,1)';
betas = paramsValues(:,2)';
gammas = paramsValues(:,3)';
lambdas = paramsValues(:,4)';

lapseFit = 'nAPLE';
gammaEQlambda = logical(false);
lapseLimits = [];
guessLimits = [];

options = [];

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions', 7)
            options = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'maxtries', 4)
            maxTries = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'rangetries', 6)
            rangeTries = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseLimits',6)
            lapseLimits = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            guessLimits = varargin{n+1};
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseFit',6)
            lapseFit = varargin{n+1};
            valid = 1;
        end                    
        if strncmpi(varargin{n}, 'gammaEQlambda',6)
            gammaEQlambda = logical(varargin{n+1});            
            valid = 1;
        end                                                
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})
        end        
    end            
end

[paramsL, LLL, exitflagL, output, funcParamsL] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsValues, PF, 'Thresholds',MC.argsAlesser,'Slopes',MC.argsBlesser,'GuessRates',MC.argsGlesser,'LapseRates',MC.argsLlesser,'searchoptions',options,'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
Try = 1;

while exitflagL == 0 && Try < maxTries    
    paramsTry = paramsValues;
    ArgsATry = MC.argsAlesser;
    ArgsBTry = MC.argsBlesser;
    ArgsGTry = MC.argsGlesser;
    ArgsLTry = MC.argsLlesser;
    if ~isstruct(rangeTries)
        multiplier = PAL_PFML_rangeTries(MC.argsAlesser,MC.argsBlesser,MC.argsGlesser,MC.argsLlesser, repmat(rangeTries,size(paramsValues,1),1));
        paramsTry = paramsValues + multiplier.*repmat(rangeTries,size(paramsValues,1),1);
    else
        if isfield(rangeTries,'rangeTries')
            multiplier = PAL_PFML_rangeTries(MC.argsAlesser,MC.argsBlesser,MC.argsGlesser,MC.argsLlesser, repmat(rangeTries.rangeTries,size(paramsValues,1),1));
            paramsTry = paramsValues + multiplier.*repmat(rangeTries.rangeTries,size(paramsValues,1),1);
        end
        if isstruct(ArgsATry) && isfield(rangeTries,'paramsValuesA');
            ArgsATry.paramsValuesA = ArgsATry.paramsValuesA + ArgsATry.paramsFreeA.*(rand(1,length(rangeTries.paramsValuesA))-.5).*rangeTries.paramsValuesA;
        end
        if isstruct(ArgsBTry) &&  isfield(rangeTries,'paramsValuesB');
            ArgsBTry.paramsValuesB = ArgsBTry.paramsValuesB + ArgsBTry.paramsFreeB.*(rand(1,length(rangeTries.paramsValuesB))-.5).*rangeTries.paramsValuesB;
        end
        if isstruct(ArgsGTry) &&  isfield(rangeTries,'paramsValuesG');
            ArgsBTry.paramsValuesG = ArgsGTry.paramsValuesG + ArgsGTry.paramsFreeG.*(rand(1,length(rangeTries.paramsValuesG))-.5).*rangeTries.paramsValuesG;
        end
        if isstruct(ArgsLTry) &&  isfield(rangeTries,'paramsValuesL');
            ArgsLTry.paramsValuesL = ArgsLTry.paramsValuesL + ArgsLTry.paramsFreeL.*(rand(1,length(rangeTries.paramsValuesL))-.5).*rangeTries.paramsValuesL;
        end
    end
    [paramsL, LLL, exitflagL, output, funcParamsL] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsTry, PF, 'Thresholds',ArgsATry,'Slopes',ArgsBTry,'GuessRates',ArgsGTry,'LapseRates',ArgsLTry,'searchoptions',options,'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    Try = Try + 1;    
end

if isstruct(MC.argsAlesser) && isstruct(MC.argsAfuller)
    temp = MC.argsAfuller;
    temp.paramsValuesA(temp.paramsFreeA == 1) = funcParamsL.paramsValuesA(temp.paramsFreeA == 1);
    MC.argsAfuller = temp;
end
if isstruct(MC.argsBlesser) && isstruct(MC.argsBfuller)
    temp = MC.argsBfuller;
    temp.paramsValuesB(temp.paramsFreeB == 1) = funcParamsL.paramsValuesB(temp.paramsFreeB == 1);
    MC.argsBfuller = temp;
end
if isstruct(MC.argsGlesser) && isstruct(MC.argsGfuller)
    temp = MC.argsGfuller;
    temp.paramsValuesG(temp.paramsFreeG == 1) = funcParamsL.paramsValuesG(temp.paramsFreeG == 1);
    MC.argsGfuller = temp;
end
if isstruct(MC.argsLlesser) && isstruct(MC.argsLfuller)
    temp = MC.argsLfuller;
    temp.paramsValuesL(temp.paramsFreeL == 1) = funcParamsL.paramsValuesL(temp.paramsFreeL == 1);
    MC.argsLfuller = temp;
end

[paramsF, LLF, exitflagF, output, funcParamsF] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsL, PF,'Thresholds',MC.argsAfuller,'Slopes',MC.argsBfuller,'GuessRates',MC.argsGfuller,'LapseRates',MC.argsLfuller,'searchoptions',options,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

Try = 1;

while exitflagF == 0 && Try < maxTries
    paramsTry = paramsL;
    ArgsATry = MC.argsAfuller;
    ArgsBTry = MC.argsBfuller;
    ArgsGTry = MC.argsGfuller;
    ArgsLTry = MC.argsLfuller;
    if ~isstruct(rangeTries)
        multiplier = PAL_PFML_rangeTries(MC.argsAfuller, MC.argsBfuller, MC.argsGfuller, MC.argsLfuller, repmat(rangeTries,size(paramsL,1),1));
        paramsTry = paramsL + multiplier.*repmat(rangeTries,size(paramsL,1),1);
    else
        if isfield(rangeTries,'rangeTries')
            multiplier = PAL_PFML_rangeTries(MC.argsAfuller, MC.argsBfuller, MC.argsGfuller, MC.argsLfuller, repmat(rangeTries.rangeTries,size(paramsL,1),1));
            paramsTry = paramsL + multiplier.*repmat(rangeTries.rangeTries,size(paramsL,1),1);
        end
        if isstruct(ArgsATry) && isfield(rangeTries,'paramsValuesA');
            ArgsATry.paramsValuesA = ArgsATry.paramsValuesA + (rand(1,length(rangeTries.paramsValuesA))-.5).*rangeTries.paramsValuesA;
        end
        if isstruct(ArgsBTry) && isfield(rangeTries,'paramsValuesB');
            ArgsBTry.paramsValuesB = ArgsBTry.paramsValuesB + (rand(1,length(rangeTries.paramsValuesB))-.5).*rangeTries.paramsValuesB;
        end
        if isstruct(ArgsGTry) && isfield(rangeTries,'paramsValuesG');
            ArgsGTry.paramsValuesG = ArgsGTry.paramsValuesG + (rand(1,length(rangeTries.paramsValuesG))-.5).*rangeTries.paramsValuesG;
        end
        if isstruct(ArgsLTry) && isfield(rangeTries,'paramsValuesL');
            ArgsLTry.paramsValuesL = ArgsLTry.paramsValuesL + (rand(1,length(rangeTries.paramsValuesL))-.5).*rangeTries.paramsValuesL;
        end
    end
    [paramsF, LLF, exitflagF, output, funcParamsF] = PAL_PFML_FitMultiple(StimLevels,NumPos, OutOfNum, paramsTry, PF,'Thresholds',ArgsATry,'Slopes',ArgsBTry,'GuessRates',ArgsGTry,'LapseRates',ArgsLTry,'searchoptions',options,'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
    Try = Try + 1;
end

TLR = 2*(LLF - LLL);
exitflag = exitflagF && exitflagL;