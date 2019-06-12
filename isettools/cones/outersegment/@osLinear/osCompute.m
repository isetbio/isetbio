function [current, interpFilters, meanCur] = osCompute(obj, cMosaic, varargin)
% Linear model computing outer seg. photocurrent from isomerizations (R*) 
%
% Syntax:
%	[current, interpFilters, meanCur] = osCompute(obj, cMosaic, varargin)
%
% Description:
%    We use  osLinear.osCompute (linear model) for experiments with a
%    uniform background, as we often find in psychophysical experiments.
%    When the images are more complex (e.g., natural scenes), use the
%    osBioPhys model, not the linear model.
%
%    When the current is computed for many trials, we can save computation
%     time by returning these values on the first call and reusing.
%
%        interpFilters - Used to save computation time in loops
%        meanCur       - Used to save computation time in loops
%
%    The linear model impulse response function is the small signal of the
%    osBioPhys model. The impulse response depends on the on mean
%    isomerization rate.
%
%    We calculate the osBioPhys current response to
%
%        * a constant abosrption rate
%        * a constant with 1 photon added to one sampling bin. 
%
%    The difference between these two signals is the photocurrent
%    impulseResponse to a photon.
%
%    The current predicted to an arbitrary stimulus is the current from the
%    mean isomerization rate plus the current from small deviations around
%    the the mean isomerization rate.
%
%        (mean current) + ...
%             convolve((absorptions - meanAbsorptions), impulseResponse) 
%
%    If the os.noiseFlag is 'random' or 'frozen', the method adds noise to
%    the current output signal. See osAddNoise() for specification and
%    validation of the noise model. See Angueyra and Rieke (2013, Nature
%    Neuroscience) for details.
%
% Inputs: 
%	 obj           - osLinear class object
%    cMosaic       - parent object of the outerSegment
%    varargin      - (Optional) Additional option variables, such as:
%           seed: (Optional) If the noise is frozen, you can send in a
%                 seed. Default is 1.
%
% Output:
%	 current       - outer segment photocurrent current in pA
%    interpFilters - Interpolated impulse response functions (to
%                    integration time samples)
%	 meanCur       - mean value of the output current
%
% Optional key/value pairs:
%    'seed',                - Random seed.
%    'interpFilters',       - The LMS linear filters to use
%    'meanCur'              - The background current
%    'absorptionsInXWFormat'- The absorptions in 2D format (space, time)
%

% History:
%    xx/xx/16  JRG/HJ/BW  ISETBIO TEAM, 2016
%    02/13/18  jnm        Formatting

%% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'outerSegment'));
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.addParameter('seed', 1, @isnumeric);
p.addParameter('interpFilters', [], @isnumeric);
p.addParameter('meanCur', [], @isnumeric);
p.addParameter('absorptionsInXWFormat', [], @isnumeric);
p.parse(obj, cMosaic, varargin{:});

% Frozen noise seed
seed = p.Results.seed;
interpFilters = p.Results.interpFilters;
meanCur = p.Results.meanCur; 
absorptionsInXWFormat = p.Results.absorptionsInXWFormat;

coneType = cMosaic.pattern;
meanRate = coneMeanIsomerizations(cMosaic, 'perSample', true, ...
    'absorptionsInXWFormat', absorptionsInXWFormat);  % R*/sample

if (~isempty(p.Results.absorptionsInXWFormat))
   tSamples = size(absorptionsInXWFormat, 2);
else
   tSamples = size(cMosaic.absorptions, 3);
end

%% Get the linear filters for the mean rate

if isempty(interpFilters) || isempty(meanCur)
    % These convert a single photon increment on mean to a photocurrent
    % impulse response function
    [lmsFilters, meanCur] = obj.linearFilters(cMosaic, ...
        'absorptionsInXWFormat', absorptionsInXWFormat, ...
        'eccentricity', obj.eccentricityDegs);
    % obj.plot('current filters', 'meancurrent', meanCur)
    
    %% Interpolate stored lmsFilters to the time base of the absorptions
    osTimeAxis = obj.timeAxis;
    absTimeAxis = cMosaic.interpFilterTimeAxis;
    interpFilters = zeros(numel(absTimeAxis), 3);
    for ii = 1:3
        % Interpolation assumes that we are accounting for the time sample
        % bin width elsewhere. Also, we extrapolate the filters with zeros
        % to make sure that they extend all the way through the absorption
        % time axis. See notes in s_matlabConv2.m for an explanation why.
        interpFilters(:, ii) = interp1(...
            osTimeAxis(:), lmsFilters(:, ii), absTimeAxis(:), 'linear', 0);
    end
end

%%  The predicted photocurrent is
if (~isempty(p.Results.absorptionsInXWFormat))
    % Already in (space, t) format
    absorptions = p.Results.absorptionsInXWFormat;
    
    % We will store the current here
    current = 0 * absorptions;
    
    % Find LMS cone indices
    nonNullConeIndices = find(cMosaic.pattern > 1);
    nonNullConeTypes = coneType(nonNullConeIndices);
    coneIndices{2} = find(nonNullConeTypes == 2);
    coneIndices{3} = find(nonNullConeTypes == 3);
    coneIndices{4} = find(nonNullConeTypes == 4);
else
    % Convert (x, y, t) to (space, t)
    [absorptions, r, c] = RGB2XWFormat(cMosaic.absorptions);  % Per sample
    % vcNewGraphWin;
    % plot(absorptions(100, :));

    % We will store the current here
    current = zeros(r*c, tSamples);
    
    % Find LMS cone indices
    coneIndices{2} = find(coneType == 2);
    coneIndices{3} = find(coneType == 3);
    coneIndices{4} = find(coneType == 4);
end

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank, so we skip
    
    % locate cones with specific type and convolve with temporal filter
    index = coneIndices{ii};
    if ~isempty(index)
        % The mean absorptions produces a mean current. This current level
        % was returned above when we calculated the lmsFilters (meanCur).
        % 
        % Here we calculate the time-varying photocurrent by  
        %   * Convolving the difference of the absorption rate from the
        %     mean absorption rate with the L, M or S filter
        %   * Adding in the mean background current
        %
        %  conv(absorptions - meanAbsorptions, lmsFilters) + meanCur
        
        % dAbsorptions is [nCones by nTime]
        dAbsorptions = absorptions(index, :) - meanRate(ii - 1);
        % The difference should be distributed around 0
        %
        %   vcNewGraphWin; hist(dAbsorptions(:));
        %   mean(dAbsorptions(:))
        
        % Convolve and  add in the mean. The general conv2 produces a new
        % time series that is longer than nTimes. We only want the
        % convolution up to the final absorption. Not really sure if we
        % want 'same' here, or we want circonv, or ... (BW).
        % tmpCurrent = conv2(dAbsorptions, ...
        %     interpFilters(:, ii - 1)', 'same') + meanCur(ii - 1);
        tmpCurrent = conv2(interpFilters(:, ii - 1)', dAbsorptions) ...
            + meanCur(ii - 1);
        % vcNewGraphWin; plot(tmpCurrent(10, :))
        
        % Store it
        current(index, :) = tmpCurrent(:, 1:tSamples);
    end
end

% Noise anyone?
switch obj.noiseFlag
    case 'none'
        % fprintf('No current noise added.\n')
    case 'random'
        % fprintf('Random noise added.\n')
        current = osAddNoise(current, 'sampTime', cMosaic.integrationTime);
    case 'frozen'
        %warning('ISETBIO:ConeMosaic:osCompute:displayFrozenNoiseSeed', ...
        %	'Frozen noise (seed %d)\n', seed);
        current = osAddNoise(current, ...
            'sampTime', cMosaic.integrationTime, 'seed', seed);
    otherwise
        error('Noise flag %s\n', obj.noiseFlag);
end

if (isempty(p.Results.absorptionsInXWFormat))
    % Reshape the current back into (x, y, t) format
    current = XW2RGBFormat(current, r, c);
    % ieMovie(current);
end

end