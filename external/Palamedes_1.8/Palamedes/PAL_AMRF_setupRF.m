%
%PAL_AMRF_setupRF  Creates structure which contains settings for and 
%   results of running fit adaptive method.
%   
%   syntax: RF = PAL_AMRF_setupRF({optional arguments})
%
%   RF = PAL_AMRF_setupRF creates and returns a structure containing
%   settings for the running fit adaptive method using default settings.
%
%   Default settings may be changed by providing pairwise arguments, where
%   the first entry of the pair is a string indicating the field to be
%   changed and the second entry provides the new value of the field.
%
%   Modifiable fields and settings (default values in curly brackets):
%   
%   'priorAlphaRange'     vector  {[-2:.01:2]}
%       Vector containing values of threshold to be considered in 
%       posterior distribution.
%
%   'prior'               vector {uniform across priorAlphaRange}
%       Prior distribution.
%
%   'beta'                positive scalar {2}
%       Slope of PF to be used in fits.
%
%   'gamma'               scalar in range [0-1] {.5}
%       Guess rate to be used in fits.
%
%   'lambda'              scalar in range [0-1] {.02}
%       Lapse rate to be used in fits.
%
%   'PF'                  inline function {@PAL_Gumbel}
%       Form of psychometric function to be assumed by Psi method.
%
%   'stopCriterion'       string {'trials'}
%       Criterion by which to terminate run. Possible values: 'trials',
%       'reversals'.
%
%   'stopRule'            positive scalar {50}
%       Value of stop criterion (trials or reversals) after which to
%       terminate run.
%
%   'xMin'                scalar or empty array {[]}
%       minimum value to be assigned to stimulus magnitude. If empty matrix
%       is assigned, range of stimulus magnitude has no lower bound.
%
%   'xMax'                scalar or empty array {[]}
%       maximum value to be assigned to stimulus magnitude. If empty matrix
%       is assigned, range of stimulus magnitude has no upper bound.
%
%   'meanmode'            string {'mean'}
%       Stimulus placement rule. Stimulus magnitude on any trial will be 
%       set to either the mean or mode of the pdf based on previous trials.
%       Possible values: 'mean', 'mode'.
%
%   Example: RF = PAL_AMRF_setupRF('stopCriterion','reversals', ...
%       'stopRule' ,12) creates a new structure using default settings for 
%       all fields except 'stopCriterion' and 'stopRule' (RF.stop will be 
%       set to 1 after 12 reversals).
%
%   In order to change settings in an existing structure, pass the existing
%   structure as the first argument. For example, given an existing
%   structure PM the call:
%       RF = PAL_AMRF_setupRF(RF, 'gamma',.25)
%   changes field gamma in the existing structure RF to .25 without
%   affecting other settings in RF.
%
%   RF's result storage fields:
%
%   'RF.xCurrent' contains stimulus magnitude to be used on current trial
%   'RF.x' stores stimulus intensities for all trials
%   'RF.response' stores responses for all trials
%   'RF.pdf' stores posterior distribution
%   'RF.mean' stores the expected value of alpha in posterior, which may be
%       used as threshold estimate
%   'RF.mode' stores the modal value of alpha in posterior, which may be
%       used as threshold estimate
%   'RF.sd' stores standard deviation of posterior, which may be used as 
%       the standard error of the threshold estimate.
%   'RF.meanUniformPrior', 'RF. modeUniformPrior', and 'RF.sdUniformPrior'
%       ignore a user-defined prior in determining mean, mode, and sd of 
%       posterior and use a rectangular prior instead.
%   'RF.stop' is used as termination flag. While stop criterion has not 
%       been reached, 'RF.stop' will equal 0, when criterion is reached, 
%       'RF.stop' will be set to 1.
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.2.0, 1.4.0, 1.6.3 (see History.m)

function RF = PAL_AMRF_setupRF(varargin)

NumOpts = length(varargin);

if mod(NumOpts,2) == 0
    RF.priorAlphaRange = -2:.01:2;
    RF.prior = ones(size(RF.priorAlphaRange));
    RF.prior = RF.prior./sum(RF.prior);
    RF.pdf = RF.prior;
    [RF.mode, RF.mean, RF.sd] = PAL_AMRF_pdfDescriptives(RF.pdf, RF.priorAlphaRange);
    RF.modeUniformPrior = [];
    RF.meanUniformPrior = [];
    RF.sdUniformPrior = [];
    RF.response = [];
    RF.stopCriterion = 'trials';
    RF.stopRule = 50;
    RF.stop = 0;
    RF.PF = @PAL_Gumbel;
    RF.beta = 2;
    RF.gamma = 0.5;
    RF.lambda = 0.02;
    RF.xMin = [];
    RF.xMax = [];
    RF.direction = [];
    RF.reversal = 0;
    RF.meanmode = 'mean';
    RF.xCurrent = RF.mean;
    RF.x = [];
    RF.xStaircase = RF.mean;    
else
    RF = varargin{1};
end
    
if NumOpts > 1
    opts(1) = cellstr('priorAlphaRange');
    opts(2) = cellstr('prior');
    opts(3) = cellstr('beta');    
    opts(4) = cellstr('gamma');
    opts(5) = cellstr('lambda');
    opts(6) = cellstr('PF');
    opts(7) = cellstr('stopCriterion');    
    opts(8) = cellstr('stopRule');
    opts(9) = cellstr('xMin');
    opts(10) = cellstr('xMax');
    opts(11) = cellstr('meanmode');
    supplied = logical(false(size(opts)));
    for opt = 1:length(opts)        
        for n = 1:2:NumOpts-mod(NumOpts,2)
            n = n+mod(NumOpts,2);
            valid = 0;
            if strncmpi(varargin{n}, opts(1),6)            
                RF.priorAlphaRange = varargin{n+1};
                valid = 1;
                supplied(1) = true;
            end
            if strcmpi(varargin{n}, opts(2))
                RF.prior = varargin{n+1};
                valid = 1;
                supplied(2) = true;
            end
            if strcmpi(varargin{n}, opts(3))
                RF.beta = varargin{n+1};
                valid = 1;
                supplied(3) = true;
            end
            if strncmpi(varargin{n}, opts(4),4)
                RF.gamma = varargin{n+1};
                valid = 1;
                supplied(4) = true;
            end
            if strncmpi(varargin{n}, opts(5),4)
                RF.lambda = varargin{n+1};
                valid = 1;
                supplied(5) = true;
            end
            if strcmpi(varargin{n}, opts(6))
                RF.PF = varargin{n+1};
                valid = 1;
                supplied(6) = true;
            end
            if strncmpi(varargin{n}, opts(7),5)
                RF.stopCriterion = varargin{n+1};
                valid = 1;
                supplied(7) = true;
            end
            if strncmpi(varargin{n}, opts(8),5)
                RF.stopRule = varargin{n+1};
                valid = 1;
                supplied(8) = true;
            end
            if strcmpi(varargin{n}, opts(9))
                RF.xMin = varargin{n+1};
                valid = 1;
                supplied(9) = true;
            end
            if strcmpi(varargin{n}, opts(10))
                RF.xMax = varargin{n+1};
                valid = 1;
                supplied(10) = true;
            end
            if strcmpi(varargin{n}, opts(11))
                RF.meanmode = varargin{n+1};
                valid = 1;
                supplied(11) = true;
            end
            if valid == 0
                warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})
            end    
        end
    end            
    if supplied(1) && ~supplied(2)
        RF.prior = ones(size(RF.priorAlphaRange));
    end
    RF.prior = RF.prior./sum(RF.prior);
    if isempty(RF.x)   %First session. Otherwise keep going with existing RF.pdf
        RF.pdf = RF.prior;
    end
    [RF.mode, RF.mean, RF.sd] = PAL_AMRF_pdfDescriptives(RF.pdf, RF.priorAlphaRange);

    if strcmpi(RF.meanmode, 'mean')
        RF.xCurrent = RF.mean;
    end
    if strcmpi(RF.meanmode, 'mode')
        RF.xCurrent = RF.mode;
    end
    if RF.xCurrent > RF.xMax
        RF.xCurrent = RF.xMax;
    end
    if RF.xCurrent < RF.xMin
        RF.xCurrent = RF.xMin;
    end
end