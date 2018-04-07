function [absorptions, current, interpFilters, meanCur] = compute(obj, oi, varargin)
%Compute the cone absorptions, possibly for multiple trials (repeats)
%
% Syntax:
%   [absorptions, current, interpFilters, meanCur] = compute(obj, oi);
%   [absorptions, current, interpFilters, meanCur] = ...
%       compute(obj, oiSequence);
%
% Description:
%    Compute the temporal Sequence of cone absorptions, which we treat as
%    isomerizations, R*. The computation can executed on
%      * a single optical image (snapshot)
%      * a single optical image with eye movements
%      * an optical image sequence (oiSequence) with eye movements.
%
%     Any additional key/value pairs will be passed on to routine the
%     computeForOISequence method, if the input argumment oi is in fact an
%     oiSequence object.
%
%     An eye movement path (emPath) can be generated using the
%     emGenSequence method, or it can be sent in as the 'emPath' variable.
%     For a single trial, the emPath is a nTime x 2 matrix of (row, col)
%     positions with respect to the cone mosaic. For the single trial case,
%     we recommend setting the coneMosaic.emPositions or using
%     emGenSequence.
%
%     When using an oiSequence input, you can execute a multiple trial
%     computation by setting the eye movement variable to a 3D array
%          emPath: (nTrials , nTime , 2)
%     In that case, we return the absorptions and possibly photocurrent for
%     nTrials in a 4D matrix:  (nTrials x nTime x row x col).
%
%     If you set 'currentFlag' to true, then the current is also returned
%     in a matrix of the same size.
%
%     The cone photon noise is computed according to the setting of the
%     obj.noiseFlag, which can be 'random', 'frozen', or 'none'. If
%     'frozen', then you can set the 'seed' parameter. The default when a
%     mosaic is created is 'random'.
%
%     The cone photocurrent is computed according to obj.os.noiseFlag,
%     which can also be set to 'random', 'frozen', or 'none', as above. The
%     default when an os object is created is 'random'.
%
% Inputs:
%    obj                 - A coneMosaic object
%    oi                  - An optical image, or oiSequence. See oiCreate
%                          for more details.
%
% Outputs:
%    absorptions         - cone photon absorptions
%    current             - cone photocurrent
%    interpFilters       - photon to photocurrent impulse response function
%    meanCur             - mean current level
%
% Optional key/value pairs:
%    'currentFlag'       - Compute photocurrent. Default false.
%    'seed'              - Seed to use when obj.noiseFlag is 'frozen'.
%                          Default 1.
%    'emPath'            - Eye movement path (Nx2 matrix) If the input is
%                          an oiSequence, the eye movement paths can be
%                          (nTrials x nTimes x 2) and multiple trials will
%                          be returned (see below for more details).
%                          Default obj.emPositions.
%    'theExpandedMosaic' - This allows you to set a larger cone mosaic to
%                          make sure that the eye movements are within the
%                          mosaic.
%
% See Also:
%    coneMosaic, computeForOISequence, emGenSequence.
%

% History:
%    xx/xx/16  JRG/BW  ISETBIO Team 2016
%    08/09/17  dhb     Working on standardizing comment format. I'm wasnt'
%                      happy with my previous pass.
%    02/26/18  jnm     Formatting

%% If an oi sequence, head that way
%
% Send to the specialized compute in that case.
if isequal(class(oi), 'oiSequence')
    [absorptions, current, interpFilters, meanCur] = ...
        obj.computeForOISequence(oi, varargin{:});
    return;
end

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('oi', @isstruct);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('seed', 1, @isnumeric);
p.addParameter('emPath', obj.emPositions, @isnumeric);
p.addParameter('theExpandedMosaic', [], @(x)(isa(x, 'coneMosaic')));
p.parse(oi, varargin{:});

currentFlag = p.Results.currentFlag;
seed = p.Results.seed;
emPath = p.Results.emPath;
theExpandedMosaic = p.Results.theExpandedMosaic;

obj.absorptions = [];
obj.current = [];

%% Set eye movement path
%
%  We do not accept multiple trials for this computational path (single
%  image). In that case you should run a loop on the compute itself.
%
if length(size(emPath)) > 2
    % If multiple trials are requested, then emPath is 
    %
    %      (nTrials x nTime x 2)
    %
    % We make a recursive call to this routine in order to compute with one
    % eye movement path at a time. The current flag is maintained.

    % Define variables and allocate space
    nTrials = size(emPath, 1);
    nTime = size(emPath, 2);
    absorptions = zeros(nTrials, obj.rows, obj.cols, nTime);
    if currentFlag, current = zeros(size(absorptions)); end

    % Do each trial
    for ii = 1:nTrials
        thisPath = squeeze(emPath(ii, :, :));
        [absorptions(ii, :, :, :), c, interpFilters, meanCur] = ...
            obj.compute(oi, 'emPaths', thisPath, ...
            'currentFlag', currentFlag);
        if currentFlag, current(ii, :, :, :) = c; end
    end

else
    
    obj.emPositions = emPath;
    
    % This code efficiently calculates the effects of eye movements. The
    % logic is this:
    %
    %   1. Make a full LMS calculation so that we know the LMS absorptions
    %   at every cone mosaic position. Only need to do this once.
    %
    %   2. Then move the eye to a position and pull out the LMS values that
    %   match the spatial pattern of the cones in the grid, but centered at
    %   the next eye movement location.
    %

    % We need a copy of the object because of eye movements.
    if (isempty(theExpandedMosaic))
        % We need an enlarged version of the cone mosaic because of the eye
        % movements. Generate it here, if the user did not provide one.
        padRows = max(abs(emPath(:, 2)));
        padCols = max(abs(emPath(:, 1)));
        theExpandedMosaic = obj.copy();
        theExpandedMosaic.pattern = zeros(obj.rows + 2 * padRows, ...
            obj.cols + 2 * padCols);
    elseif isa(theExpandedMosaic, 'coneMosaic')
        % OK, we are passed theExpandedMosaic.
        % Set the current path and integrationTime and use it.
        theExpandedMosaic.emPositions = obj.emPositions;
        theExpandedMosaic.integrationTime = obj.integrationTime;
        theExpandedMosaic.absorptions = [];
        padRows = round((theExpandedMosaic.rows - obj.rows) / 2);
        padCols = round((theExpandedMosaic.cols - obj.cols) / 2);
    end

    % Compute full LMS noise free absorptions
    absorptions = theExpandedMosaic.computeSingleFrame(oi, ...
        'fullLMS', true);

    % Deal with eye movements
    absorptions = obj.applyEMPath(absorptions, 'emPath', emPath, ...
        'padRows', padRows, 'padCols', padCols);

    % Add photon noise to the whole volume
    switch obj.noiseFlag
        case {'frozen', 'random'}
            if (isa(obj, 'coneMosaicHex'))
                % Only call photonNoise on non-null cones for a hex mosaic.
                nonNullConeIndices = find(obj.pattern > 1);
                timeSamples = size(absorptions, 3);
                absorptions = reshape(permute(absorptions, [3 1 2]), ...
                    [timeSamples size(obj.pattern, 1) * ...
                    size(obj.pattern, 2)]);
                absorptionsCopy = absorptions;
                absorptions = absorptions(:, nonNullConeIndices);

                % Add noise
                absorptionsCopy(:, nonNullConeIndices) = ...
                    obj.photonNoise(absorptions, ...
                    'noiseFlag', obj.noiseFlag, 'seed', seed);
                absorptions = permute(reshape(absorptionsCopy, ...
                    [timeSamples size(obj.pattern, 1), ...
                    size(obj.pattern, 2)]), [2 3 1]);
                clear 'absorptionsCopy'
            else % Rectangular mosaic
                % Add noise
                absorptions = obj.photonNoise(absorptions, ...
                    'noiseFlag', obj.noiseFlag, 'seed', seed);
            end
        case {'none'}
            % No noise
        otherwise
            error('Invalid noise flag passed');
    end

    % Set the absorptions in the object.
    obj.absorptions = absorptions;

    %% Compute photocurrent if requested
    current = [];
    interpFilters = [];
    meanCur = [];
    if currentFlag
        if size(obj.absorptions, 3) == 1
            disp(['Absorptions are a single frame. No current to ' ...
                'calculate.'])
            return;
        else
            [current, interpFilters, meanCur] = obj.os.osCompute(obj);
            obj.current = current;
        end
    end
end

end
