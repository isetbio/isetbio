%
%PAL_AMPM_setupPM  Creates structure which contains settings for and 
%   results of Kontsevich & Tyler's (1999) psi adaptive method and 
%   variations on it (see Prins (2013) or 
%   www.palamedestoolbox.org/psimarginal.html)
%   
%   syntax: PM = PAL_AMPM_setupPM({optional arguments})
%
%   PM = PAL_AMPM_setupPM creates and returns a structure containing 
%       settings for the psi method adaptive method using default settings.
%
%   Default settings may be changed by providing pairwise arguments, where
%   the first entry of the pair is a string indicating the field to be
%   changed and the second entry provides the new value of the field.
%   Modifiable fields and settings (default values in curly brackets):
%   
%   'priorAlphaRange'     vector  {[-2:.05:2]}
%       Vector containing values of threshold to be considered in posterior
%       distribution.
%   
%   'priorBetaRange'      vector  {[-1:.05:1]}
%       Vector containing log (base 10) transformed values of slope to be 
%       considered in posterior distribution.
%
%   'priorGammaRange'     scalar or vector  {[.5]}
%       Vector containing values of threshold to be considered in posterior
%       distribution.
%
%   'priorLambdaRange'     scalar or vector  {[.02]}
%       Vector containing values of threshold to be considered in posterior
%       distribution.
%      
%   'stimRange'           vector  {[-1:.1:1]}
%       Possible stimulus values to be considered on each trial.
%
%   'prior'               matrix  {uniform}
%       prior should have size: length(priorAlphaRange) x
%       length(priorBetaRange) x length(priorGammaRange) x
%       length(priorLambdaRange)
%       But: Trailing singleton dimensions are ignored, e.g., if
%       length(priorGammaRange) and length(priorLambdaRange) are both 1,
%       prior is of size: length(priorAlphaRange) x length(priorBetaRange).
%       User can specify a prior distribution using this option. For
%       example of use, see PAL_AMPM_Demo.
%
%   'PF'                  inline function {@PAL_Gumbel}
%       Form of psychometric function to be assumed by Psi method.
%
%   'numTrials'           positive integer {50}
%       Length of run in terms of number of trials
%
%   'gammaEQlambda'       logical {false}
%       When gammaEQlambda is set to true, gamma and lambda are both
%       assumed to result from lapses. Any entry for 'priorGammaRange' will
%       be ignored. Since it is assumed that gamma = lambda, a single 
%       parameter will be estimated for both, the prior used will be that 
%       set for 'priorLambdaRange'.
%
%   'marginalize'         vector {[]} (or string)
%       Allows users to marginalize parameters out of the posterior
%       distribution before the to-be-minimized expected entropy is
%       calculated. This allows one to include, say, the lapse rate in the
%       posterior distribution without the method having the explicit goal 
%       of minimizing entropy with respect to the lapse rate. Roughly
%       speaking, marginalized parameters will be estimated only insofar as
%       doing so is the most optimal manner in which to reduce entropy in 
%       the non-marginalized parameters in the posterior distribution.
%       To marginalize parameters, pass the string 'threshold', 'slope', 
%       'guess', or 'lapse' after the string 'marginalize'. To marginalize
%       more than one parameter, repeatedly call PAIR of arguments (e.g., 
%       'marginalize', 'lapse','marginalize','guess'). Parameters to be 
%       marginalized may also be specified by passing a vector 
%       containing a numerical code for each of the parameters to be added 
%       to marginalize list (1: threshold, 2: slope, 3: guess rate, 4: 
%       lapse rate). In order to remove a parameter from the list of to be 
%       marginalized parameters, use the negative symbol (i.e., -). For 
%       example, the call PM = PAL_AMPM_setupPM(PM,'marginalize',
%       '-threshold') removes the threshold from PM's marginalize list (if 
%       it was indeed on there). PM = PAL_AMPM_setupPM(PM,'marginalize',
%       [-1]) does the same. Example 2: The call: PM = 
%       PAL_AMPM_setupPM(PM,'marginalize',[-1 2 3 -4]) will result in the 
%       slope and guess rate to be marginalized, while the threshold and 
%       lapse rate will not be marginalized, whatever the previous 
%       marginalize list was. In order to delete all previous entries, pass 
%       empty vector. Note that parameters may be added to or removed from
%       list at any point during a trial run. For more information see: 
%       www.palamedestoolbox.org/psimarginal.html
%
%   Example: PM = PAL_AMPM_setupPM('numTrials',100, 'PriorLambdaRange', ...
%   [0:.01:.1],'marginalize','lapse') creates a new structure using default 
%   settings for all fields except numTrials, priorLambdaRange, and 
%   marginalize. Setting 'PriorLambdaRange' to [0:.01:.1] will mean the 
%   method will include the lapse rate as a third dimension in the
%   posterior distribution. Marginalizing the lapse rate means that on each
%   trial the Psi method will place the stimulus such as to minimize 
%   expected entropy in the posterior across threshold and slope with the 
%   lapse rate marginalized. 
%
%
%   In order to change settings in an existing structure, pass the existing
%   structure as the first argument. For example, given an existing
%   structure 'PM' the call:
%       PM = PAL_AMPM_setupPM(PM, 'gamma',.25)
%   changes the field gamma in the existing structure 'PM' to .25 without
%   affecting other settings in existing structure 'PM.'.
%
%   PM's result storage fields:
%
%   'PM.x' stores stimulus intensities for all trials
%   'PM.response' stores responses for all trials (positive (correct,
%       'greater than'): 1, negative: 0)
%   'PM.pdf' stores posterior distribution
%   'PM.threshold' stores threshold estimates after each trial (marginal 
%       expected value of alpha in posterior)
%   'PM.slope' stores log slope estimates after each trial (marginal 
%       expected value of log beta in posterior)
%   'PM.guess' stores guess rate estimates after each trial (marginal 
%       expected value of gamma in posterior)
%   'PM.lapse' stores lapse rate estimates after each trial (marginal 
%       expected value of lambda in posterior)
%   'PM.seThreshold' stores standard error of threshold estimate (marginal 
%       standard deviation of alpha in posterior).
%   'PM.seSlope' stores standard error of log slope estimate (marginal 
%       standard deviation of log beta in posterior).
%   'PM.seGuess' stores standard error of guess rate estimate (marginal 
%       standard deviation of log beta in posterior).
%   'PM.seLapse' stores standard error of lapse rate estimate (marginal 
%       standard deviation of log beta in posterior).
%   'PM.thresholdUniformPrior', 'PM.slopeUniformPrior', 
%       'PM.guessUniformPrior', 'PM.lapseUniformPrior', 
%       'PM.seThresholdUniformPrior', 'PM.seSlopeUniformPrior',
%       'PM.seGuessUniformPrior', and 'PM.seLapseUniformPrior' ignore 
%       user-defined prior and determine estimates using a uniform prior 
%       instead.
%   'PM.stop' is used as termination flag. While stop criterion has not 
%       been reached, 'PM.stop' will equal 0, when criterion is reached, 
%       'PM.stop' will be set to 1.
%
% References:
%
% Kontsevich, L.L. & Tyler, C.W. (1999). Bayesian adaptive estimation of 
%   psychometric slope and threshold. Vision Research, 39, 2729–2737.
%
% Prins, N. (2013). The psi-marginal adaptive method: how to give nuisance 
%   parameters the attention they deserve (no more, no less). Journal of
%   Vision, 13(7):3, 1-17. doi: 10.1167/13.7.3 
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.1.1, 1.2.0, 1.4.0, 1.5.0, 1.6.0, 1.6.1, 
%   1.6.3 (see History.m)

function PM = PAL_AMPM_setupPM(varargin)

NumOpts = length(varargin);

message = 'As of Palamedes 1.5.0, alpha values vary across rows and beta values vary ';
message = [message 'across columns of a 2-D (alpha x beta) PM.pdf (in earlier versions this was '];
message = [message 'vice versa). This will affect behavior of code written for previous versions '];
message = [message 'ONLY if a user-supplied prior was used and/or if code performed operations '];
message = [message 'on PM.pdf (e.g., if PM.pdf was plotted using for example ''surf''). Code written for '];
message = [message 'earlier versions of Palamedes can be made compatible simply '];
message = [message 'by transposing any user-supplied prior (i.e., prior = prior'';) before passing it to '];
message = [message 'PAL_AMPM_setupPM and transposing PM.pdf before plotting it '];
message = [message '(e.g., surf(PM.pdf'') instead of surf(PM.pdf)). Visit '];
message = [message 'www.palamedestoolbox.org/pal_ampm_incompatibility.html.'];
warning('PALAMEDES:AMPM_setupPM:priorTranspose',message);
warning('off','PALAMEDES:AMPM_setupPM:priorTranspose');

if mod(NumOpts,2) == 0

    PM.priorAlphaRange = -2:.05:2;
    PM.priorBetaRange = -1:.05:1;
    PM.priorGammaRange = 0.5;
    PM.priorLambdaRange = 0.02;
    PM.gammaEQlambda = logical(false);
    PM.stimRange = -1:.1:1;
    [PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(PM.priorAlphaRange,PM.priorBetaRange,PM.priorGammaRange,PM.priorLambdaRange);
    PM.PF = @PAL_Gumbel;
    PM.LUT = PAL_AMPM_CreateLUT(PM.priorAlphaRange, PM.priorBetaRange, PM.priorGammaRange, PM.priorLambdaRange, PM.stimRange, PM.PF,PM.gammaEQlambda);
    PM.prior = ones(size(PM.priorAlphas));
    PM.prior = PM.prior./sum(sum(sum(sum(PM.prior))));
    PM.pdf = PM.prior;
    [PM.posteriorTplus1givenSuccess, PM.posteriorTplus1givenFailure, pSuccessGivenx] = PAL_AMPM_PosteriorTplus1(PM.pdf, PM.LUT); 
    ExpectedEntropy = PAL_Entropy(PM.posteriorTplus1givenSuccess,4).*pSuccessGivenx + PAL_Entropy(PM.posteriorTplus1givenFailure,4).*(1-pSuccessGivenx);
    [minEntropy, PM.I] = min(squeeze(ExpectedEntropy));
    PM.xCurrent = PM.stimRange(PM.I);
    PM.x = PM.xCurrent;
    PM.numTrials = 50;
    PM.response = [];
    PM.stop = 0;
    PM.marginalize = [];
else 
    PM = varargin{1};
end

PM.firstsession = length(PM.x) == 1;

if NumOpts > 1
    opts(1) = cellstr('priorAlphaRange');
    opts(2) = cellstr('priorBetaRange');
    opts(3) = cellstr('priorGammaRange');
    opts(4) = cellstr('priorLambdaRange');
    opts(5) = cellstr('stimRange');
    opts(6) = cellstr('prior');
    opts(7) = cellstr('PF');
    opts(8) = cellstr('numTrials');
    opts(9) = cellstr('gamma');             %for compatibility with older usage
    opts(10) = cellstr('lambda');           %for compatibility with older usage
    opts(11) = cellstr('gammaEQlambda');
    opts(12) = cellstr('marginalize');
    supplied = logical(false(size(opts)));

    for n = 1:2:NumOpts-mod(NumOpts,2)
        n = n+mod(NumOpts,2);
        valid = 0;
        if strncmpi(varargin{n}, opts(1),6)            
            PM.priorAlphaRange = varargin{n+1};                
            valid = 1;
            supplied(1) = true;
        end
        if strncmpi(varargin{n}, opts(2),6)            
            PM.priorBetaRange = varargin{n+1};
            valid = 1;
            supplied(2) = true;
        end
        if strncmpi(varargin{n}, opts(3),6) || (strncmpi(varargin{n}, opts(9),5) && ~strncmpi(varargin{n}, opts(11), 6))
            PM.priorGammaRange = varargin{n+1};                
            valid = 1;
            supplied(3) = true;
        end
        if strncmpi(varargin{n}, opts(4),6) || strncmpi(varargin{n}, opts(10),5)
            PM.priorLambdaRange = varargin{n+1};
            valid = 1;
            supplied(4) = true;
        end
        if strncmpi(varargin{n}, opts(5),4)            
            PM.stimRange = varargin{n+1};
            valid = 1;
            supplied(5) = true;
        end
        if strcmpi(varargin{n}, opts(6))
            PM.prior = varargin{n+1};
            PM.prior = PM.prior./sum(sum(sum(sum(PM.prior))));
            PM.pdf = PM.prior;
            valid = 1;
            supplied(6) = true;
        end
        if strcmpi(varargin{n}, opts(7))
            PM.PF = varargin{n+1};
            valid = 1;
            supplied(7) = true;
        end            
        if strncmpi(varargin{n}, opts(8),4)
            PM.numTrials = varargin{n+1};
            valid = 1;
            supplied(8) = true;
        end
        if strncmpi(varargin{n}, opts(11),6)
            PM.gammaEQlambda = logical(varargin{n+1});
            valid = 1;
            supplied(11) = true;
        end
        if strncmpi(varargin{n}, opts(12),6)
            switch PAL_whatIs(varargin{n+1})
                case 0
                    PM.marginalize = [];
                case 1
                    for index = 1:length(varargin{n+1})
                        if varargin{n+1}(index) > 0 && isempty(find(PM.marginalize == varargin{n+1}(index)))
                            PM.marginalize = [PM.marginalize varargin{n+1}(index)];
                        end
                        if varargin{n+1}(index) < 0
                            PM.marginalize = PM.marginalize(PM.marginalize ~= abs(varargin{n+1}(index)));
                        end
                    end
                case 2
                    if strncmpi(varargin{n+1},'threshold',5) && isempty(find(PM.marginalize == 1))
                        PM.marginalize = [PM.marginalize 1];
                    end
                    if strncmpi(varargin{n+1},'slope',5) && isempty(find(PM.marginalize == 2))
                        PM.marginalize = [PM.marginalize 2];
                    end
                    if strncmpi(varargin{n+1},'guess',5) && isempty(find(PM.marginalize == 3))
                        PM.marginalize = [PM.marginalize 3];
                    end
                    if strncmpi(varargin{n+1},'lapse',5) && isempty(find(PM.marginalize == 4))
                        PM.marginalize = [PM.marginalize 4];                    
                    end                
                    if strncmpi(varargin{n+1},'-threshold',5)
                        PM.marginalize = PM.marginalize(PM.marginalize ~= 1);
                    end
                    if strncmpi(varargin{n+1},'-slope',5)
                        PM.marginalize = PM.marginalize(PM.marginalize ~= 2);
                    end
                    if strncmpi(varargin{n+1},'-guess',5)
                        PM.marginalize = PM.marginalize(PM.marginalize ~= 3);
                    end
                    if strncmpi(varargin{n+1},'-lapse',5)
                        PM.marginalize = PM.marginalize(PM.marginalize ~= 4);                    
                    end                
            end
            valid = 1;
            supplied(12) = true;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end        
    end
    if PM.gammaEQlambda == true;
        PM.priorGammaRange = 0; %value will be ignored
    end
    if supplied(1) || supplied(2) || supplied(3) || supplied(4)
        [PM.priorAlphas, PM.priorBetas, PM.priorGammas, PM.priorLambdas] = ndgrid(PM.priorAlphaRange,PM.priorBetaRange,PM.priorGammaRange,PM.priorLambdaRange);
        if ~supplied(6)
            PM.prior = ones(size(PM.priorAlphas));
            PM.prior = PM.prior./sum(sum(sum(sum(PM.prior))));
            if PM.firstsession == 1    %First session. Otherwise keep going with existing PM.pdf
                PM.pdf = PM.prior;
            end
        end
    end
    
    if supplied(1) || supplied(2) || supplied(3) || supplied(4) || supplied(5) || supplied(7) || supplied(12)
        PM.LUT = PAL_AMPM_CreateLUT(PM.priorAlphaRange,PM.priorBetaRange,PM.priorGammaRange,PM.priorLambdaRange,PM.stimRange,PM.PF,PM.gammaEQlambda);

        [PM, expectedEntropy] = PAL_AMPM_expectedEntropy(PM);
        
        [minEntropy, PM.I] = min(squeeze(expectedEntropy));
        PM.xCurrent = PM.stimRange(PM.I);
        PM.x(length(PM.x)) = PM.xCurrent;
    end
    if PM.firstsession == 1
        PM.x(1) = PM.xCurrent;
    end
end

if PM.firstsession == 1
    PM.threshold = []; 
    PM.slope = [];
    PM.guess = [];
    PM.lapse = [];
    PM.seThreshold = [];
    PM.seSlope = [];
    PM.seGuess = [];
    PM.seLapse = [];
    
    PM.thresholdUniformPrior = [];
    PM.slopeUniformPrior = [];
    PM.guessUniformPrior = [];    
    PM.lapseUniformPrior = [];
    PM.seThresholdUniformPrior = [];
    PM.seSlopeUniformPrior = [];
    PM.seGuessUniformPrior = [];
    PM.seLapseUniformPrior = [];
end