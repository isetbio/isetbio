%
%PAL_AMRF_updateRF  Updates structure which contains settings for and 
%   results of running fit adaptive method.
%
%   syntax: RF = PAL_AMRF_updateRF(RF, amplitude, response)
%
%   After having created a structure 'RF' using PAL_AMRF_setupRF, use 
%   something akin to the following loop to control stimulus intensity
%   during experimental run:
%
%   while ~RF.stop
%       
%       %Present trial here at stimulus magnitude in 'RF.xCurrent'
%       %and collect response (1: correct/greater than, 0: incorrect/
%       %smaller than)
%
%       RF = PAL_AMRF_updateRF(RF, amplitude, response); %update RF 
%                                   %structure based on stimulus magnitude 
%                                   %and response                                    
%    
%   end
%
%
%Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.4.0, 1.6.3 (see History.m)


function RF = PAL_AMRF_updateRF(RF, amplitude, response)

trial = length(RF.response)+1;
RF.x(trial) = amplitude;
RF.response(trial) = response;

if trial == 1
    RF.xStaircase(trial) = RF.x(trial);
    if response == 1
        RF.direction = -1;        
    else
        RF.direction = 1;
    end
end

RF.pdf = PAL_AMRF_pdfUpdate(RF.pdf, RF.priorAlphaRange, RF.beta, RF.gamma, RF.lambda, RF.x(trial), response, RF.PF);
[RF.mode, RF.mean, RF.sd] = PAL_AMRF_pdfDescriptives(RF.pdf, RF.priorAlphaRange);
if strcmpi(RF.meanmode,'mean')
    if (RF.mean > RF.xCurrent && RF.direction == 1) || (RF.mean < RF.xCurrent && RF.direction == -1)
        RF.reversal(trial) = 0;
    end
    if RF.mean > RF.xCurrent && RF.direction == -1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = 1;
    end
    if RF.mean < RF.xCurrent && RF.direction == 1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = -1;
    end
    RF.xCurrent = RF.mean;
end
if strcmpi(RF.meanmode,'mode')
    if (RF.mode > RF.xCurrent && RF.direction == 1) || (RF.mode < RF.xCurrent && RF.direction == -1)
        RF.reversal(trial) = 0;
    end
    if RF.mode > RF.xCurrent && RF.direction == -1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = 1;
    end
    if RF.mode < RF.xCurrent && RF.direction == 1
        RF.reversal(trial) = sum(RF.reversal~=0) + 1;
        RF.direction = -1;
    end
    RF.xCurrent = RF.mode;
end
RF.xStaircase(trial+1) = RF.xCurrent;

if (strncmpi(RF.stopCriterion,'reversals',4) && sum(RF.reversal~=0) == RF.stopRule)||(strncmpi(RF.stopCriterion,'trials',4) && trial == RF.stopRule)
    RF.stop = 1;
    [RF.modeUniformPrior, RF.meanUniformPrior, RF.sdUniformPrior] = PAL_AMRF_pdfDescriptives(RF.pdf./RF.prior, RF.priorAlphaRange);
end