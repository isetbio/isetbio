function current = osCompute(obj, pRate, coneType, varargin)
% Compute the response of the outer segments using the linear model 
%
%    current = osCompute(obj, pRate, coneType, varargin)
%
% We use the linear model for cases in which there is a uniform background,
% as we typically find in psychophysical experiments.  When the images are
% more complex (e.g., natural scenes) you should use the osBioPhys model,
% not the linear model.
%
% The impulse response function of the linear model calculated here matches
% the linear osBioPhys model given that the observer is in a steady-state
% on the mean background (specified by the three cone pRate values).
%
% The basic idea is this.  When we have a mean field with a specific photon
% absorption rate, we get the base current returned from the osBioPhys()
% model.  We then get the temporal impulse response function that we expect
% when we add 1 photon to the mean level.  The predicted current, then, is
%
%   (absorptions - meanAbsorptions) * impulseResponse + (base current)
% 
% This converts isomerizations (R*) to outer segment current (pA). If the
% noiseFlag is set to true, this method adds noise to the current output
% signal. See Angueyra and Rieke (2013, Nature Neuroscience) for details.
%
% TODO:
%   pRate -> meanRates (one for each cone type)
%     These can be calculated here from the obj.absorptions.
%     So the pRate can go away.
%
% Inputs: 
%   obj      - osLinear class object
%   pRate    - photon absorption rate in R*/sec
%   coneType - cone type matrix, 1 for blank, 2-4 for LMS respectively
%
% Optional input (key-value pairs)
%   append   - logical, compute new or to append to existing. When append
%              is true, the computed current is appended to the existing
%              photocurrent data in the object. The returned current value
%              only contains photocurrent for input pRate.
% 
% Outputs:
%   current  - outer segment current in pA
% 
% JRG/HJ/BW, ISETBIO TEAM, 2016

% We don't think anyone sends in a sensor any more.
%
% check pRate type for backward compatibility
% if isstruct(pRate) && isfield(pRate, 'type') && ...
%         strcmp(pRate.type, 'sensor')
%     warning('The input is a sensor, should update to use coneMosaic.');
%     obj.osSet('timestep', sensorGet(pRate, 'time interval'))
%     if notDefined('coneType')
%         current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
%             sensorGet(pRate, 'cone type'));
%     else
%         current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
%             sensorGet(pRate, 'cone type'), coneType, varargin{:});
%     end
%     in the old code, we return obj instead of current
%     current = obj.osSet('cone current signal', current);
%     return
% end

% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osLinear'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);  % Comes from coneMosaic parent

% To remove and write a script to check model compatibility with Juan's stuff
p.addParameter('linearized', true, @islogical);

p.parse(obj, pRate, coneType, varargin{:});
linearized = p.Results.linearized;
coneType   = p.Results.coneType;

nHistFrames = 0;

lConeIndices = find(coneType == 2);
mConeIndices = find(coneType == 3);
sConeIndices = find(coneType == 4);

%% We only need a 3-vector pRate for three cone types.
%
% Thats because we are assuming a uniform field as in a psychophysical or
% physiological experiment for the linear case.


pRateXW = RGB2XWFormat(pRate);
lConeAbsorptions = pRateXW(lConeIndices,:);
mConeAbsorptions = pRateXW(mConeIndices,:);
sConeAbsorptions = pRateXW(sConeIndices,:);

if ~isempty(lConeAbsorptions); 
    lConeMean = mean(lConeAbsorptions(:)); 
else 
    lConeMean = 0;
end;
if ~isempty(mConeAbsorptions); 
    mConeMean = mean(mConeAbsorptions(:)); 
else
    mConeMean = 0;
end;
if ~isempty(sConeAbsorptions); 
    sConeMean = mean(sConeAbsorptions(:)); 
else
    sConeMean = 0;
end;
pMean = [lConeMean mConeMean sConeMean];

% This will get moved to a specific test of the other filters
%    lmsFilters = obj.generateLinearFilters(mean(pMean(:))); % linear filters

%% This is the place where we get the linear filters given the mean rate

% scaleFactor = 0.11143; % from physiology experiments, see coneIRFtutorial.m
% 
% Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)
% Ib = pMean; %[7131 6017 1973]; % R* per sec due to background adapting field (one for each cone, L, M, S)
% % adjust this to specific experiment
% % Ib = [2250 2250 2250];
% gain_dark = 0.32;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
% gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field

% call this bioPhysLinearFilters() to produce the impulse response
% functions for the specific mean rates
lmsFilter = obj.generateBioPhysFilters('meanRate', pMean, varargin{:}); % linear filters

% lmsFilter0 = lmsFilter(1);
% lmsFilter = scaleFactor*(lmsFilter - lmsFilter0);
% 
% % scale IRF to reflect amplitude at chosen background
% % using Weber adaptation equation above and gainRatio derived from it
% newGain = gainRatio .* gain_dark;
% oldGain = max(lmsFilter);
% IRFScaleFactor = newGain ./ oldGain;
% 
% lmsFilters = (IRFScaleFactor'*lmsFilter' - ones(3,size(lmsFilter,1))*lmsFilter0)';

%% These convert a single photon increment on mean to a photocurrent IRF

obj.lmsConeFilter = lmsFilters;

%% We need to ask whether we can't put this into the units

% This converts the mean absorption rate into the base current
maxCur = 0.01*20.5^3; % Dunn & Rieke (2007, Nature)
meanCur = maxCur * (1 - 1./(1 + 45000./pMean));

%%  The predicted photocurrent is
%
%  conv(absorptions - meanAbsorptions,lmsFilters) + baseCurrent
%


% First entry is trial.  We are showing only the first trial here.
[absorptions, r, c] = RGB2XWFormat(pRate);

current = zeros(r*c, size(pRate, 3));

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    % pull out the linear filter for current cone type.
    filter = lmsFilters(:, ii-1);
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        curData = conv2(absorptions(index, :), filter') - meanCur(ii-1);
        current(index, :) = curData(:, nHistFrames+1+(1:size(pRate, 3)));
    end
end

% record only recent history in obj
newStart = max(size(pRate, 3) - length(filter) + 1, 1);
pRate = pRate(:, :, newStart:end);

% reshape the output signal matrix.
current = XW2RGBFormat(current, r, c);

% Add noise
% The osAddNoise function expects and input to be isomerization rate.
% This is handled properly because the params has the time sampling
% rate included.
if osGet(obj,'noiseFlag') == 1
    disp('Current noise added.')
    params.sampTime = obj.timeStep;
    current = osAddNoise(current, params);
else
    disp('No current noise added.')
end


obj.coneCurrentSignal = current;

end