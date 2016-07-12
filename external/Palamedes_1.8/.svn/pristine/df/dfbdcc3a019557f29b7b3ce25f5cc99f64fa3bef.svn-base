%
%PAL_AMUD_updateUD  Updates structure which contains settings for and 
%   results of up/down adaptive method.
%
%   syntax: UD = PAL_AMUD_updateUD(UD, response)
%
%   After having created a structure 'UD' using PAL_AMUD_setupUD, use 
%   something akin to the following loop to control stimulus intensity
%   during experimental run:
%
%   while ~UD.stop
%       
%       %Present trial here at stimulus magnitude in 'UD.xCurrent'
%       %and collect response (1: correct/greater than, 0: incorrect/
%       %smaller than)
%
%       UD = PAL_AMUD_updateUD(UD, response); %update UD structure based 
%                                             %on response                                    
%    
%   end
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.4.0, 1.4.5 (see History.m)

function UD = PAL_AMUD_updateUD(UD, response)

trial = length(UD.x);
UD.response(trial) = response;

if trial == 1
    UD.xStaircase(trial) = UD.x(trial);
    if response == 1
        UD.direction = -1;        
    else
        UD.direction = 1;
    end
end
    
if response == 1
    UD.d = UD.d + 1;
    if UD.d == UD.down || max(UD.reversal) < 1
        UD.xStaircase(trial+1) = UD.xStaircase(trial)-UD.stepSizeDown;
        if UD.xStaircase(trial+1) < UD.xMin && strcmp(UD.truncate,'yes')
            UD.xStaircase(trial+1) = UD.xMin;
        end
        UD.u = 0;
        UD.d = 0;
        UD.reversal(trial) = 0;
        if UD.direction == 1
            UD.reversal(trial) = sum(UD.reversal~=0) + 1;
        else
            UD.reversal(trial) = 0;
        end
        UD.direction = -1;
    else
        UD.xStaircase(trial+1) = UD.xStaircase(trial);
    end    
else
    UD.u = UD.u + 1;
    if UD.u == UD.up || max(UD.reversal) < 1
        UD.xStaircase(trial+1) = UD.xStaircase(trial)+UD.stepSizeUp;
        if UD.xStaircase(trial+1) > UD.xMax && strcmp(UD.truncate,'yes')
            UD.xStaircase(trial+1) = UD.xMax;
        end
        UD.u = 0;
        UD.d = 0;
        UD.reversal(trial) = 0;
        if UD.direction == -1
            UD.reversal(trial) = sum(UD.reversal~=0) + 1;
        else
            UD.reversal(trial) = 0;
        end
        UD.direction = 1;
    else
        UD.xStaircase(trial+1) = UD.xStaircase(trial);
    end    
end    

if strncmpi(UD.stopCriterion,'reversals',4) && sum(UD.reversal~=0) == UD.stopRule
    UD.stop = 1;
end
if strncmpi(UD.stopCriterion,'trials',4) && trial == UD.stopRule
    UD.stop = 1;
end
if ~UD.stop
    UD.x(trial+1) = UD.xStaircase(trial+1);
    if UD.x(trial+1) > UD.xMax
        UD.x(trial+1) = UD.xMax;
    elseif UD.x(trial+1) < UD.xMin
        UD.x(trial+1) = UD.xMin;
    end
    UD.xCurrent = UD.x(trial+1);
end