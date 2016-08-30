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
    obj.osSet('timestep', sensorGet(pRate, 'time interval'));
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
p.addParameter('append', false, @islogical);

p.parse(obj, pRate, coneType, varargin{:});
isAppend = p.Results.append;

% init parameters
if ~isAppend, obj.absHistory = []; end % clean up stored state
if isempty(obj.absHistory)
    nHistFrames = 0;
    obj.absHistory = pRate;
else
    nHistFrames = size(obj.absHistory, 3);
    obj.absHistory = cat(3, obj.absHistory, pRate);
end

% generate temporal filters
pMean = mean(obj.absHistory(:));
lmsFilters = obj.generateLinearFilters(pMean); % linear filters

maxCur = 0.01*20.5^3; % Angueyra & Rieke (2013, Nature)
meanCur = maxCur * (1 - 1/(1 + 45000/pMean));

[absorptions, r, c] = RGB2XWFormat(obj.absHistory);
current = zeros(r*c, size(pRate, 3));
pRateRS = RGB2XWFormat(pRate);
% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    % pull out the linear filter for current cone type.
    filter = lmsFilters(:, ii-1);
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        curData = conv2(absorptions(index, :), filter') - meanCur;
        current(index, :) = curData(:, nHistFrames+1+(1:size(pRate, 3)));
        
%         if size(pRateRS,2) > size(filter,1)
%             filterZP = [repmat(filter',[size(pRateRS(index,:),1) 1]) zeros(size(pRateRS(index,:),1),-size(filter,1)+size(pRateRS,2))];
%             pRateZP = pRateRS(index,:);
%         else
%             filterZP = filter;
%             pRateZP = [pRateRS zeros(size(filter,1)-size(pRateRS,2),size(pRateRS,1))];
%         end
%         curData = ifft(fft(filterZP').*fft(pRateZP'))';        
%         current(index, :) = curData;
    end
end

% record only recent history in obj
newStart = max(size(obj.absHistory, 3) - length(filter) + 1, 1);
obj.absHistory = obj.absHistory(:, :, newStart:end);

% reshape the output signal matrix.
current = XW2RGBFormat(current, r, c);

% Add noise
% The osAddNoise function expects and input to be isomerization rate.
% This is handled properly because the params has the time sampling
% rate included.
if osGet(obj,'noiseFlag') == 1
    params.sampTime = obj.timeStep;
    current = osAddNoise(current, params);
end
if isAppend
    obj.coneCurrentSignal = cat(3, obj.coneCurrentSignal, current);
else
    obj.coneCurrentSignal = current;
end

end