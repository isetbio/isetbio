function response = compute(obj, varargin)
% Spatially and temporally filter the cone current to create bipolar output
%
% Syntax:
%    The bipolar object handles multiple trials, and these are returned
%    separately for the center and surround of the biplar model when use
%    have output arguments, as below.
%
%   [~, bpNTrialsCenter, bpNTrialsSurround] = ...
%           bp.compute(cMosaic, 'nTrialsInput', cMosaicCurrentNTrials);
%
%   @bipolarMosaic.compute(varargin);
%
% Description:
%    The bipolars act as a spatial-temporal function that converts the cone
%    photocurrent into bipolar current that is delivered to the retinal
%    ganglion cells. The cone mosaic is the input, and it contains a
%    photocurrent time series (coneMosaic.current).
%
%    The bipolar response is computed with a separable spatio-temporal
%    calculation (first spatial filtering, then temporal filtering) of the
%    cone current to create the bipolar output.
%
%    The critical bipolar mosaic parameters are set up and stored in the
%    @bipolarMosaic object itself. This routine computes the responses
%    based on those settings. See the class definition for the list of
%    parameters.
%
%    The principal decision is whether the bipolar transformation is linear
%    or includes a rectification. This is controlled by the obj.rectifyType
%    parameter. The rectification happens at the end of the function, after
%    both spatial and temporal filtering.
%
% Input:
%    obj          - A bipolar cone mosaic object
%
% Optional Key/Value Pairs:
%    nTrialsInput - The number of trials
%
% Output:
%    Response     - The response of the bipolar. N.B. The center and
%                   surround responses are calculated and stored in the
%                   mosaic as responseCenter and responseSurround
%
% References:
%   *Meister option* -
%       The spatial filtering is followed by a temporal filter that is
%       selected in order to match the impulse response that expect to find
%       at the RGC input.
%
%   *Chichilnisky option* -
%       The cMosaicCurrentNTrials has dimensions (nTrials, row, col,
%       nFrames), as returned by cmosaic.computeCurrent. The full bipolar
%       response simply add bpNTrialsCenter + bpNTrialsSurround.
%
% Notes:
%    * Only the last trial reponse is stored in the object itself.
%    * This feature has not been tested or much used - on our list
%    * TODO: This could be transformed to be like the rgc computation where
%      each cell has its own spatial RF. Not there yet, just using
%      convolutions.
%    * TODO: We should set this up as a structure that we use to implement
%      the anatomical rules. Let's send in a struct that defines the
%      anatomical rules (e.g., anatomyRules) with slots that implement the
%      kind of stuff listed here.
%    * Anatomical rules:
%      - off Diffuse, on Diffuse, and on Midget - Receive no S cone input
%      - offMidget - keep S cones but scale connection strength down 75%
%      - onSBC - S cone inputs to center, only L/M cone inputs to surround
%
% See Also:
%    Bipolar.m. Wiki <https://github.com/isetbio/isetbio/wiki/bipolar>
%

%% History:
% 5/2016 JRG (c) isetbio team
%
%    10/19/17  jnm  Comments & Formatting

% Examples:
%{
   % ETTBSkip - Example is broken. Variables undefined. Remove this line when fixed.
   [~, bpNTrialsCenter, bpNTrialsSurround] = ...
           bp.compute(cMosaic, 'nTrialsInput', cMosaicCurrentNTrials);
%}

%% parse input parameters
p = inputParser;
p.addRequired('obj', @(x) (isa(x, 'bipolarMosaic')));
addParameter(p, 'coneTrials',  [], @isnumeric);
%%%
% parse - no options at this opint
p.parse(obj, varargin{:});

coneTrials = p.Results.coneTrials;
cmosaic    = obj.input;

if isempty(cmosaic.current)
    error('No cone photocurrent. Use cmosaic.computeCurrent.');
end

if ~isempty(coneTrials),  nTrials = size(coneTrials, 1);
else,                     nTrials = 1;
end

%% Spatial filtering and subsampling
% If the input includes multiple trials, we run across all the trials here.
for iTrial = 1:nTrials
    % This places the cone 3D matrix into a coneNumber x time matrix
    if ~isempty(coneTrials)
        osSig = RGB2XWFormat(squeeze(coneTrials(iTrial, :, :, :)));
    else
        osSig = RGB2XWFormat(cmosaic.current);
    end

    %% Enforce anatomical rules on cone to bipolar connections
    switch obj.cellType
        case{'offdiffuse', 'ondiffuse', 'onmidget'}
            osSigCenter   = osSig;
            osSigSurround = osSig;

            % Remove S cone input for these types of bipolars
            %
            % Find the locations indices of the different cone types
            [~, ~, S] = coneTypeLocations(cmosaic, 'format', 'index');

            % Zero the photocurrent of the S cones. Do this for both the
            % center and the surround.
            minval = min(osSig(:));
            osSigCenter(S(:), :)   = minval*ones(size(osSigCenter(S, :)));
            osSigSurround(S(:), :) = minval*ones(size(osSigCenter(S, :)));

        case{'offmidget'}
            % Keep S cone input for off Midget but only weight by 0.25
            % Find the locations (row, col) of the different cone types
            [~, ~, S] = coneTypeLocations(cmosaic, 'format', 'index');
            minval = min(osSig(:));

            osSigCenter   = osSig;
            osSigCenter(S, :)   = 0.25*(osSigCenter(S, :)-minval)+minval;

            osSigSurround   = osSig;
            osSigSurround(S, :) = 0.25*(osSigSurround(S, :)-minval)+minval;

        case{'onsbc'}
            % Set L and M cones to zero in SBC center, set S cones to zero
            % in SBC surround.
            % Find the indices of the different cone types
            [L, M, S] = coneTypeLocations(cmosaic, 'format', 'index');

            % This is one long vector of L, M cone indices
            LM = [L; M];

            % Find the effectively zero outer segment signal for this
            % mosaic
            minval = min(osSig(:));

            % When the center is an LM cone, make all of the time steps in
            % the center the smallest value
            osSigCenter       = osSig;
            osSigCenter(LM, :) = minval*ones(size(osSigCenter(LM, :)));

            % Put effectively zero S-cone signals into the surround
            osSigSurround      = osSig;
            osSigSurround(S, :) = minval*ones(size(osSigSurround(S, :)));

        otherwise
            error('Unrecognized bipolar mosaic type %s\n', obj.cellType);
    end

    % Put the data back into RGB format, like RGB2XW()
    sz = size(cmosaic.current);
    osSigCenter   = XW2RGBFormat(osSigCenter, sz(1), sz(2));
    osSigSurround = XW2RGBFormat(osSigSurround, sz(1), sz(2));
    % cmosaic.window;
    % vcNewGraphWin; ieMovie(osSigCenter);

    %% Spatial filtering
    % Full spatial convolution for every frame. The kernel is only 2D
    % which is why we have a space-only convolution.
    bipolarCenter   = ieSpaceTimeFilter(osSigCenter, obj.sRFcenter);
    bipolarSurround = ieSpaceTimeFilter(osSigSurround, obj.sRFsurround);
    % vcNewGraphWin; ieMovie(bipolarCenter);
    % vcNewGraphWin; ieMovie(bipolarSurround);

    % Pull out the samples at the cell locations. It works here because
    % they are evenly spaced (stride). If we have jitter, we need another
    % approach.
    strideRow = abs(obj.cellLocation(1, 2, 1) - obj.cellLocation(1, 1, 1));
    strideCol = abs(obj.cellLocation(2, 1, 2) - obj.cellLocation(1, 1, 2));

    bipolarCenter   = bipolarCenter(1:strideRow:end, 1:strideCol:end, :);
    bipolarSurround = bipolarSurround(1:strideRow:end, 1:strideCol:end, :);

    %% Temporal filtering
    % Reshape the data for the temporal convolution
    [bipolarCenter, row, col] = RGB2XWFormat(bipolarCenter);
    bipolarSurround = RGB2XWFormat(bipolarSurround);

    % This is the impulse response filter
    bipolarFilt = bipolarFilter(obj, cmosaic);

    % If we wanted to rectify the signal, we could do it here
    % obj.rectify(input, 'rType', {hw, fw, none})
    % obj.responseCenter   = ...
    %    obj.rectificationCenter(bipolarOutputLinearCenter);
    % obj.responseSurround = ...
    %    obj.rectificationSurround(bipolarOutputLinearSurround);
    %

    %% Rectification and temporal convolution issues
    % Rectification - not tested or analyzed
    %
    % We have in the past shifted the bipolar response levels to a minimum
    % of zero. That is arbitrary and produces higher contrast signals. It
    % might be OK because we have no real units on the bipolar current.
    % Or, maybe we should leave them alone. Anyway, here we are shifting
    % them to a min of zero.

    % bipolarCenter = obj.rectificationCenter(bipolarCenter ...
    %     - (min(bipolarCenter')'*ones(1, size(bipolarCenter, 2))));
    bipolarCenter = bipolarCenter - (min(bipolarCenter, [], 2) ...
        * ones(1, size(bipolarCenter, 2)));
    tmpCenter = conv2(bipolarFilt, bipolarCenter);
    % vcNewGraphWin; tmp = XW2RGBFormat(tmpCenter, row, col); ieMovie(tmp);

    % Rectification
    % Not fully tested or analyzed -
    % bipolarSurround = obj.rectificationSurround(bipolarSurround ...
    %     - (min(bipolarSurround')'*ones(1, size(bipolarSurround, 2))));
    bipolarSurround = bipolarSurround-(min(bipolarSurround, [], 2) ...
        * ones(1, size(bipolarSurround, 2)));
    tmpSurround = conv2(bipolarFilt, bipolarSurround);
    % vcNewGraphWin;
    % tmp = XW2RGBFormat(tmpSurround, row, col); ieMovie(tmp);

    if ~isempty(coneTrials)
        if iTrial == 1
            nTrialsCenter = zeros([nTrials, size(XW2RGBFormat(...
                tmpCenter(:, 1:cmosaic.tSamples), row, col))]);
            nTrialsSurround = zeros([nTrials, size(XW2RGBFormat(...
                tmpSurround(:, 1:cmosaic.tSamples), row, col))]);
        end

        nTrialsCenter(iTrial, :, :, :) = XW2RGBFormat(...
            tmpCenter(:, 1:cmosaic.tSamples), row, col);
        nTrialsSurround(iTrial, :, :, :) = XW2RGBFormat(...
            tmpSurround(:, 1:cmosaic.tSamples), row, col);
    else
        nTrialsCenter   = 0;
        nTrialsSurround = 0;
    end

    if iTrial == nTrials
        % Store the last trial in the object
        obj.responseCenter   = XW2RGBFormat(tmpCenter(:, ...
            1:size(cmosaic.current, 3)), row, col);
        obj.responseSurround = XW2RGBFormat(tmpSurround(:, ...
            1:size(cmosaic.current, 3)), row, col);
    end
end

response = nTrialsCenter - nTrialsSurround;
end
