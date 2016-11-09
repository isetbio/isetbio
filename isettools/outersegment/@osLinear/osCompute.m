function current = osCompute(obj, pRate, coneType, varargin)
% Compute the linear filter response of the outer segments. 
%
%    current = osCompute(obj, pRate, coneType, varargin)
%
% This converts isomerizations (R*) to outer segment current (pA). If the
% noiseFlag is set to true, this method adds noise to the current output
% signal. See Angueyra and Rieke (2013, Nature Neuroscience) for details.
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

% check pRate type for backward compatibility
if isstruct(pRate) && isfield(pRate, 'type') && ...
        strcmp(pRate.type, 'sensor')
    warning('The input is a sensor, should update to use coneMosaic.');
    obj.osSet('timestep', sensorGet(pRate, 'time interval'))
    if notDefined('coneType')
        current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
            sensorGet(pRate, 'cone type'));
    else
        current = obj.osCompute(sensorGet(pRate, 'photon rate'), ...
            sensorGet(pRate, 'cone type'), coneType, varargin{:});
    end
    % in the old code, we return obj instead of current
    current = obj.osSet('cone current signal', current);
    return
end

% parse inputs
p = inputParser; p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osLinear'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);
p.addParameter('linearized', true, @islogical);

p.parse(obj, pRate, coneType, varargin{:});
linearized = p.Results.linearized;

nHistFrames = 0;

lConeIndices = find(coneType == 2);
mConeIndices = find(coneType == 3);
sConeIndices = find(coneType == 4);

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

if ~linearized
    lmsFilters = obj.generateLinearFilters(mean(pMean(:))); % linear filters
else
    scaleFactor = 0.11143; % from physiology experiments, see coneIRFtutorial.m
    
    Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)
    Ib = pMean; %[7131 6017 1973]; % R* per sec due to background adapting field (one for each cone, L, M, S)
    % adjust this to specific experiment
    % Ib = [2250 2250 2250];
    gain_dark = 0.32;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
    gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field
    
    lmsFilter = obj.generateBioPhysFilters('meanRate', pMean, varargin{:}); % linear filters
    lmsFilter0 = lmsFilter(1);
    lmsFilter = scaleFactor*(lmsFilter - lmsFilter0);
    
    % scale IRF to reflect amplitude at chosen background
    % using Weber adaptation equation above and gainRatio derived from it
    newGain = gainRatio .* gain_dark;
    oldGain = max(lmsFilter);
    IRFScaleFactor = newGain ./ oldGain;
    
    lmsFilters = (IRFScaleFactor'*lmsFilter' - ones(3,size(lmsFilter,1))*lmsFilter0)';
end

obj.lmsConeFilter = lmsFilters;
maxCur = 0.01*20.5^3; % Angueyra & Rieke (2013, Nature)
meanCur = maxCur * (1 - 1./(1 + 45000./pMean));

% First entry is trial.  We are showing only the first trial here.
[absorptions, r, c] = RGB2XWFormat(pRate);

% pRateRS = RGB2XWFormat(pRate);
% if size(pRateRS,2) > size(lmsFilters(:, 1),1)
current = zeros(r*c, size(pRate, 3));
% else    
%     current = zeros(r*c, size(lmsFilters(:, 1), 1));
% end

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    % pull out the linear filter for current cone type.
    filter = lmsFilters(:, ii-1);
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        curData = conv2(absorptions(index, :), filter') - meanCur(ii-1);
        current(index, :) = curData(:, nHistFrames+1+(1:size(pRate, 3)));
        
%         if size(pRateRS,2) > size(filter,1)
%             filterZP = [repmat(filter',[size(pRateRS(index,:),1) 1]) zeros(size(pRateRS(index,:),1),-size(filter,1)+size(pRateRS,2))];
%             pRateZP = pRateRS(index,:);
%         else
%             filterZP = repmat(filter',[size(pRateRS(index,:),1) 1]);
%             pRateZP = [pRateRS(index,:) zeros(size(pRateRS(index,:),1),size(filter,1)-size(pRateRS(index,:),2))];
%         end
%         curData = ifft(fft(filterZP').*fft(pRateZP'))';        
%         current(index, :) = curData;
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