%
%PAL_PFML_FitMultiple   Fit psychometric functions to several conditions 
%   simultaneously using a Maximum Likelihood criterion.
%
%Syntax: [paramsValues LL exitflag output funcParams numParams] = ...
%           PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%           paramsValues, PF,{optional arguments})
%
%Input: 
%   'StimLevels': matrix containing stimulus levels used. Each row of the
%       matrix corresponds to a condition.
%
%   'NumPos': matrix containing for each of the entries of 'StimLevels' the 
%       number of trials on which a positive response (e.g., 'yes' or 
%       'correct') was given.
%
%   'OutOfNum': matrix containing for each of the entries of 'StimLevels' 
%       the total number of trials. (May contain zeros in case not all 
%       conditions utilized identical numbers of stimulus levels).
%
%   'paramsValues': Matrix containing parametervalues. Each row should 
%       correspond to a condition. Matrix should have four columns 
%       [threshold slope guess-rate lapse-rate]. Entries should indicate 
%       initial guesses for free parameters, fixed values for fixed 
%       parameters. 'paramsValues' may also be a 1 x 4 vector, in which
%       case the entries will be used as initial guesses for all
%       conditions.
%
%   'PF': psychometric function to be fitted. Passed as an inline function.
%       Options include:    
%           @PAL_Logistic
%           @PAL_Weibull
%           @PAL_Gumbel (i.e., log-Weibull)
%           @PAL_Quick
%           @PAL_logQuick
%           @PAL_CumulativeNormal
%           @PAL_Gumbel
%           @PAL_HyperbolicSecant
%
%Output:
%   'paramsValues': Matrix containing values of fitted and fixed 
%       parameters of the psychometric functions.
%
%   'LL': Log likelihood associated with the collective fit.
%
%   'exitflag': 1 indicates a succesful fit, 0 indicates fit did not
%       converge (trying again using new initial guesses might help).
%
%   'output': message containing some information concerning fitting
%       process.
%
%   'funcParams': structure containing parameter estimates of custom-
%       defined reparametrizations (if any). Type 'help
%       PAL_PFML_CustomDefine' for more information.
%
%   'numParams': number of free parameters.
%
%Unless specified otherwise, PAL_PFML_FitMultiple will fit all
%   conditions separately treating threshold and slope as free parameters
%   and the guess-rate and lapse-rate as fixed parameters. However, other
%   fitting schemes may be specified by user. For all of the four
%   parameters of a psychometric function the following options may be
%   specified (example of use is given below):
%       'fixed': makes the parameter fixed for all conditions. Parameter
%           will be set to value specified by user in 'paramsValues'.
%       'constrained': parameter is free, but constrained to have identical
%           value across all conditions.
%       'unconstrained': parameter is free, and may take on different 
%           values across the different conditions.
%       User may also provide a 'contrast matrix'. A contrast
%           matrix has as many columns as there are conditions and may have
%           up to as many rows as there are conditions. Given n conditions,
%           a mxn contrast matrix transforms the n parameters into m
%           new parameters, each a weighted sum of the original
%           parameters. PAL_PFML_FitMultiple will find best fitting values
%           for the m transformed parameters then transform these m back
%           into the original metric. Example: let's say there are four
%           conditions and thus also four thresholds. The contrast matrix
%
%           [ 1  1  1  1;
%            -1 -1  1  1]
%
%           will transform the four thresholds {a1 thru a4} into two new
%           parameters: b1 = (1)a1 + (1)a2 + (1)a3 + (1)a4 and b2 = (-1)a1
%           + (-1)a2 + (1)a3 + (1)a4. Best-fitting values of b1 and b2 will
%           be found and transformed back into 'a' values. This scheme will
%           allow conditions 1 and 2 on the one hand and conditions 3 and 4
%           on the other to have different thresholds. However, the
%           thresholds will be constrained to be identical between 
%           condition 1 and 2 as well as between condition 3 and 4. See
%           example below. Judd, McClelland and Ryan [1] gives a
%           gentle introduction to the use of contrasts to define models.
%       Finally, user may custom-define constraints on parameter values. 
%           Type 'help PAL_PFML_CustomDefine' for more information.
%
%Different fitting schemes may be specified using the optional argument
%   'lapseFits'. Default value is 'nAPLE' (non-Asymptotic Performance Lapse 
%   Estimation), in which all free parameters are estimated jointly in the 
%   manner outlined by e.g., Wichmann & Hill (2001). In case the highest 
%   stimulus value in 'StimLevels' is so high that it can be assumed that 
%   errors at this intensity can only be due to lapses, alternative schemes 
%   may be specified. 'jAPLE' (joint APLE) assumes that errors at
%   highest 'StimLevel' are due exclusively to lapses, and observations at 
%   this intensity contribute only to estimate of lapse rate. Observations 
%   at other values in 'StimLevels' contribute to all free parameters 
%   (including lapse rate). 'iAPLE' (isolated APLE) is identical to 
%   'jAPLE', except that the lapse rate estimate is based on observations 
%   made at the highest 'StimLevel' only. Other parameters are estimated 
%   from observations made at the other values of 'StimLevel' using lapse 
%   rate estimate obtained at highest 'StimLevel' as fixed value. Options 
%   may be used in conjunction with 'constrained' and 'uncnstrained' 
%   options on 'lapserates' (see above), as well as with contrast-defined 
%   constraints on lapse rates. See example below.
%
%The guess rate and lapse rate parameter can be constrained to be equal, as
%   would be appropriate, for example, in a bistable percept task. To
%   accomplish this, use optional argument 'gammaEQlambda', followed by a 
%   1. Both the guess rate and lapse rate parameters will be fit according
%   to options set for the lapse rate parameter. Value for guess rate in
%   'paramsValues' will be ignored.
%
%User may constrain the lapse rate to fall within limited range using 
%   the optional argument 'lapseLimits', followed by two-element vector 
%   containing lower and upper limit respectively. See full example below.
%
%User may constrain the guess rate to fall within limited range using 
%   the optional argument 'guessLimits', followed by two-element vector 
%   containing lower and upper limit respectively.
%
%PAL_PFML_FitMultiple uses Nelder-Mead Simplex method to find the maximum 
%   in the likelihood function. The default search options may be changed 
%   by using the optional argument 'SearchOptions' followed by an options 
%   structure created using options = PAL_minimize('options'). See example 
%   of usage below. For more information type PAL_minimize('options','help').
%
%
%Full examples:
%
%   options = PAL_minimize('options');  %decrease tolerance (i.e., increase
%   options.TolX = 1e-09;              %precision). This is a good idea,
%   options.TolFun = 1e-09;            %especially in high-dimension
%   options.MaxIter = 10000;           %parameter space.
%   options.MaxFunEvals = 10000;
%
%   StimLevels = [-2:1:2; -2:1:2; -2:1:2; -2:1:2];
%   NumPos = [158 177 222 268 275;
%               196 232 282 351 380;
%               159 161 218 253 283;
%               213 245 296 356 375];
%   OutOfNum = [300 300 300 300 300;
%                400 400 400 400 400;
%                300 300 300 300 300;
%                400 400 400 400 400];
%   PF = @PAL_Logistic;
%   paramsValues = [0 1.5 .5 .04];
%
% %Fitted using defaults:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF,'searchOptions',options) 
%
%returns:
%   
%   paramsFitted =
% 
%    -0.0706    1.3951    0.5000    0.0400
%     0.0955    1.7277    0.5000    0.0400
%     0.1631    1.5478    0.5000    0.0400
%    -0.1318    1.4193    0.5000    0.0400
%
% %Fitted assuming that errors made at the highest stimulus intensity can
% only be due to lapses and that lapse rates are identical between first
% and second condition and also identical between third and fourth 
% conditions:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF,'lapserates',[1 1 1 1; 1 1 -1 -1],'lapseFit',...
%       'iAPLE','searchOptions',options)  
%
%returns:
%
%   paramsFitted =
% 
%    -0.1714    1.7485    0.5000    0.0643
%     0.0132    1.8032    0.5000    0.0643
%     0.0952    1.5537    0.5000    0.0600
%    -0.2119    1.5285    0.5000    0.0600
%
%One more example:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
%       paramsValues, PF, 'thresholds', [1 1 1 1; -1 -1 1 1],'slopes',...
%       'fixed','guessrates','fixed','lapserates','constrained',...
%       'lapselimits',[0 .05],'searchOptions',options)
%
%returns:
%
%   paramsFitted =
% 
%     0.0172    1.5000    0.5000    0.0393
%     0.0172    1.5000    0.5000    0.0393
%     0.0092    1.5000    0.5000    0.0393
%     0.0092    1.5000    0.5000    0.0393
%
%[1] Judd, C.M., McClelland, G.H., & Ryan, C.S. (2008). Data analysis: A 
%   model comparison approach. Routledge.
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.0.2, 1.1.0, 1.3.0, 1.4.0, 1.4.2, 1.6.0,
%   1.6.3 (see History.m)

function [paramsValues, LL, exitflag, output, funcParams, numParams] = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, paramsValues, PF, varargin)

options = [];
lapseLimits = [];
guessLimits = [];
lapseFit = 'nAPLE';
gammaEQlambda = logical(false);
funcParams = PAL_PFML_setupParametrizationStruct;

if size(paramsValues,1) == 1
    paramsValues = repmat(paramsValues,[size(StimLevels,1) 1]);
end

FM.argsA = 'unconstrained';
FM.argsB = 'unconstrained';
FM.argsG = 'fixed';
FM.argsL = 'fixed';

if ~isempty(varargin)
    NumOpts = length(varargin);
    n = 1;
    while n <= NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'SearchOptions',7)
            options = varargin{n+1};
            valid = 1;
            add = 2;
        end
        if strncmpi(varargin{n}, 'Thresholds',5)
            FM.argsA = varargin{n+1};
            add = 2;
            valid = 1;
        end
        if strcmpi(varargin{n}, 'Slopes')
            FM.argsB = varargin{n+1};
            add = 2;
            valid = 1;
        end
        if strncmpi(varargin{n}, 'GuessRates',6)
            FM.argsG = varargin{n+1};
            add = 2;
            valid = 1;
        end
        if strncmpi(varargin{n}, 'LapseRates',6)
            FM.argsL = varargin{n+1};
            add = 2;
            valid = 1;
        end
        if strncmpi(varargin{n}, 'lapseLimits',6)
            lapseLimits = varargin{n+1};
            valid = 1;
            add = 2;
        end
        if strncmpi(varargin{n}, 'guessLimits',6)
            guessLimits = varargin{n+1};
            valid = 1;
            add = 2;
        end
        if strncmpi(varargin{n}, 'lapseFit',6)
            lapseFit = varargin{n+1};
            valid = 1;
            add = 2;
        end        
        if strncmpi(varargin{n}, 'gammaEQlambda',6)
            gammaEQlambda = logical(varargin{n+1});            
            valid = 1;
            add = 2;
        end                        
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end        
        n = n + add;
    end            
end

if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && ((PAL_whatIs(FM.argsL) == 2 && strncmpi(FM.argsL,'fix',3)) || isempty(FM.argsL))
    warning('PALAMEDES:invalidOption','Lapses rates are not free: ''LapseLimits'' argument ignored');
    lapseFit = 'nap';
end
if (strncmpi(lapseFit,'iap',3) || strncmpi(lapseFit,'jap',3)) && PAL_whatIs(FM.argsL) == 4
    warning('PALAMEDES:invalidOption','Lapse rates custom-reparameterized: ''lapseFit'' argument ignored');
    lapseFit = 'nap';
end

if gammaEQlambda
    FM.argsG = [];
    if ~isempty(guessLimits)
        warning('PALAMEDES:invalidOption','Guess rates constrained to equal lapse rates: ''guessLimits'' ignored');
        guessLimits = [];
    end
end

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos,OutOfNum);

if ~PAL_PFML_IndependentFit(FM)

    [thetas, thetasID, FM] = PAL_PFML_PtoT(paramsValues, FM);
    numParams = length(thetas);

    if strncmpi(lapseFit,'iap',3)
        for cond=1:size(NumPos,1);
            Index = find(StimLevels(cond,:) == max(StimLevels(cond,OutOfNum(cond,:)~=0)));
            NonLapsesN(cond)=NumPos(cond,Index);
            totalN(cond)=OutOfNum(cond,Index);
        end
        if gammaEQlambda
            for cond=1:size(NumPos,1);
                Index = find(StimLevels(cond,:) == min(StimLevels(cond,OutOfNum(cond,:)~=0)));
                NonLapsesN(cond)=NonLapsesN(cond)+(OutOfNum(cond,Index)-NumPos(cond,Index));
                totalN(cond)=totalN(cond) + OutOfNum(cond,Index);
            end            
        end
        switch PAL_whatIs(FM.argsL)            
            case 1
                paramsValues(:,4) = 1-((NonLapsesN/FM.argsL)*FM.argsL)./((totalN/FM.argsL)*FM.argsL);
            case 2                
                if strncmpi(FM.argsL,'unc',3)
                    paramsValues(:,4) = 1-(NonLapsesN./totalN);
                else
                    paramsValues(:,4) = sum(NonLapsesN)./sum(totalN);
                end
                if ~isempty(lapseLimits)
                    paramsValues(:,4) = max(paramsValues(:,4),lapseLimits(1));
                    paramsValues(:,4) = min(paramsValues(:,4),lapseLimits(2));
                end
        end        
        FM.argsL = [];
        thetas = thetas(thetasID~=4);
        thetasID = thetasID(thetasID~=4);        
    end    


    [thetas, negLL, exitflag, output] = PAL_minimize(@PAL_PFML_negLLMultiple, thetas, options, thetasID, paramsValues, StimLevels, NumPos, OutOfNum, FM, PF, 'lapseLimits', lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);

    paramsValues = PAL_PFML_TtoP(paramsValues, thetas, thetasID, FM);
    
    if isstruct(FM.argsA)
        funcParams.funcA = FM.argsA.funcA;
        funcParams.paramsValuesA = ones(1,length(FM.argsA.paramsValuesA));
        funcParams.paramsValuesA(FM.argsA.paramsFreeA == 1) = thetas(thetasID == 1);
        funcParams.paramsValuesA(FM.argsA.paramsFreeA == 0) = FM.argsA.paramsValuesA(FM.argsA.paramsFreeA == 0);
        funcParams.paramsFreeA = FM.argsA.paramsFreeA;        
    end
    if isstruct(FM.argsB)
        funcParams.funcB = FM.argsB.funcB;
        funcParams.paramsValuesB = ones(1,length(FM.argsB.paramsValuesB));
        funcParams.paramsValuesB(FM.argsB.paramsFreeB == 1) = thetas(thetasID == 2);
        funcParams.paramsValuesB(FM.argsB.paramsFreeB == 0) = FM.argsB.paramsValuesB(FM.argsB.paramsFreeB == 0);
        funcParams.paramsFreeB = FM.argsB.paramsFreeB;        
    end
    if isstruct(FM.argsG)
        funcParams.funcG = FM.argsG.funcG;
        funcParams.paramsValuesG = ones(1,length(FM.argsG.paramsValuesG));
        funcParams.paramsValuesG(FM.argsG.paramsFreeG == 1) = thetas(thetasID == 3);
        funcParams.paramsValuesG(FM.argsG.paramsFreeG == 0) = FM.argsG.paramsValuesG(FM.argsG.paramsFreeG == 0);
        funcParams.paramsFreeG = FM.argsG.paramsFreeG;        
    end
    if isstruct(FM.argsL)
        funcParams.funcL = FM.argsL.funcL;
        funcParams.paramsValuesL = ones(1,length(FM.argsL.paramsValuesL));
        funcParams.paramsValuesL(FM.argsL.paramsFreeL == 1) = thetas(thetasID == 4);
        funcParams.paramsValuesL(FM.argsL.paramsFreeL == 0) = FM.argsL.paramsValuesL(FM.argsL.paramsFreeL == 0);
        funcParams.paramsFreeL = FM.argsL.paramsFreeL;        
    end
    
    LL = -1*negLL;

else
    freeA = ~((ischar(FM.argsA) && strncmpi(FM.argsA,'fix',3)) || (isnumeric(FM.argsA) && isempty(FM.argsA)));
    freeB = ~((ischar(FM.argsB) && strncmpi(FM.argsB,'fix',3)) || (isnumeric(FM.argsB) && isempty(FM.argsB)));
    freeG = ~((ischar(FM.argsG) && strncmpi(FM.argsG,'fix',3)) || (isnumeric(FM.argsG) && isempty(FM.argsG)));
    freeL = ~((ischar(FM.argsL) && strncmpi(FM.argsL,'fix',3)) || (isnumeric(FM.argsL) && isempty(FM.argsL)));
    paramsFree = [freeA freeB freeG freeL];
    numParams = sum(paramsFree).*size(StimLevels,1);
    for cond = 1:size(StimLevels,1)        
        paramsFreeVals = paramsValues(cond,paramsFree == 1);
        paramsFixedVals = paramsValues(cond,paramsFree == 0);
        
        if isempty(paramsFreeVals)
            negLL(cond) = PAL_PFML_negLL(paramsFreeVals, paramsFixedVals, paramsFree, StimLevels(cond,:), NumPos(cond,:), OutOfNum(cond,:), PF,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
            exitflagCond(cond) = 1;
            output = ['No parameters estimated: None were free to vary'];
        else            
            if strncmpi(lapseFit,'iap',3)
                 len = length(NumPos(cond,OutOfNum(cond,:)~=0));
                if ~gammaEQlambda
                    lambda = 1 - NumPos(cond,len)./OutOfNum(cond,len);
                    if ~isempty(lapseLimits)
                        lambda = min(lambda,lapseLimits(2));
                        lambda = max(lambda,lapseLimits(1));
                    end                    
                    paramsFixedVals(length(paramsFixedVals)+1) = lambda; %lapse rate
                    paramsFreeVals = paramsFreeVals(1:length(paramsFreeVals)-1);    %set lapse rate estimate as fixed value
                    paramsFree(4) = 0;
                    [paramsFreeVals, negLL(cond), exitflagCond(cond), output] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels(cond,1:len-1), NumPos(cond,1:len-1), OutOfNum(cond,1:len-1), PF,'gammaEQlambda',gammaEQlambda,'guessLimits', guessLimits);        
                    negLL(cond) = negLL(cond) - log((1 - lambda).^NumPos(cond,len)) - log(lambda.^(OutOfNum(cond,len)-NumPos(cond,len)));
                else
                    lambda = 1 - (NumPos(cond,len)+(OutOfNum(cond,1)-NumPos(cond,1)))./(OutOfNum(cond,len)+OutOfNum(cond,1));
                    if ~isempty(lapseLimits)
                        lambda = min(lambda,lapseLimits(2));
                        lambda = max(lambda,lapseLimits(1));
                    end                    
                    paramsFixedVals(length(paramsFixedVals)+1) = lambda; %lapse rate
                    paramsFixedVals(length(paramsFixedVals)) = lambda; %guess rate
                    paramsFreeVals = paramsFreeVals(1:length(paramsFreeVals)-1);    %set lapse rate estimate as fixed value
                    paramsFree(4) = 0;
                    [paramsFreeVals, negLL(cond), exitflagCond(cond), output] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels(cond,2:len-1), NumPos(cond,2:len-1), OutOfNum(cond,2:len-1), PF,'gammaEQlambda',gammaEQlambda,'lapseLimits', lapseLimits);        
                    negLL(cond) = negLL(cond) - log((1 - lambda).^(NumPos(cond,len)+(OutOfNum(cond,1)-NumPos(cond,1)))) - log(lambda.^(OutOfNum(cond,len)-NumPos(cond,len)+NumPos(cond,1)));
                end            
                 
            else
                [paramsFreeVals, negLL(cond), exitflagCond(cond), output(cond)] = PAL_minimize(@PAL_PFML_negLL, paramsFreeVals, options, paramsFixedVals, paramsFree, StimLevels(cond,:), NumPos(cond,:), OutOfNum(cond,:), PF, 'lapseLimits',lapseLimits,'guessLimits', guessLimits,'lapseFit',lapseFit,'gammaEQlambda',gammaEQlambda);
            end
        end
        paramsValues(cond,paramsFree == 1) = paramsFreeVals;
        paramsValues(cond,paramsFree == 0) = paramsFixedVals;
        paramsFree(4) = freeL;
    end
    exitflag = (sum(exitflagCond)==size(StimLevels,1));
    LL = -sum(negLL);
end
if gammaEQlambda
    paramsValues(:,3) = paramsValues(:,4);
end
