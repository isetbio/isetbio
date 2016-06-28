function current = osCompute(obj, pRate, coneType, varargin)
% Compute the linear filter response of the outer segments. 
%
%    obj = osCompute(obj, sensor, varargin)
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
% Outputs:
%   current  - outer segment current in pA
% 
% JRG/HJ/BW, ISETBIO TEAM, 2016

% parse inputs
p = inputParser; p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osLinear'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);

p.parse(obj, pRate, coneType, varargin{:});

% init parameters
lmsFilters = obj.generateLinearFilters(mean(pRate(:))); % linear filters
nFrames = size(pRate, 3);

maxCur = 0.01*20.5^3; % Angueyra & Rieke (2013, Nature)
meanCur = maxCur * (1 - 1/(1 + 45000/mean(pRate(:))));

[pRate, r, c] = RGB2XWFormat(pRate);
current = zeros(size(pRate));

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    % pull out the linear filter for current cone type.
    filter = lmsFilters(:, ii-1);
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        curData = conv2(pRate(index, :), filter') - meanCur;
        current(index, :) = curData(:, 2:nFrames+1);
    end
end

% % Reshape the output signal matrix.
current = XW2RGBFormat(current, r, c);

% Add noise
% The osAddNoise function expects and input to be isomerization rate.
% This is handled properly because the params has the time sampling
% rate included.
if osGet(obj,'noiseFlag') == 1
    params.sampTime = obj.timeStep;
    current = osAddNoise(current, params);
end
obj.osSet('cone current signal', current);

end