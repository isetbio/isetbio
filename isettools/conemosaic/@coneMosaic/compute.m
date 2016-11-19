function [absorptions, current, varargout] = compute(obj, oi, varargin)
% Compute the pattern of cone absorptions and possibly the photocurrent
%
%  [absorptions, current] = cMosaic.compute(oi);
%
% Inputs:
%   oi  - optical image, or oiSequence.  See oiCreate for more details
%
% Optional inputs:
%   currentFlag  - logical, whether to compute photocurrent
%   newNoise     - logical, whether to use new random seed
%   emPath       - eye movement path in nx2 matrix. This
%                  parameter shadows obj.emPositions and is
%                  required when append is true
%
% Outputs:
%   absorptions  - cone photon absorptions
%   current      - cone photocurrent
%
% HJ ISETBIO Team 2016

%% If an oi sequence, head that way

% Send to the specialized compute in that case.
if isequal(class(oi),'oiSequence')
    % Deleted various returns that pertainin to current time axis.  In the
    % modern age, the current and absorption time axes are the same and
    % derived from the integrationTime and number of samples in the
    % coneMosaic object. (BW)
    [absorptions, current] = obj.computeForOISequence(oi,varargin{:});
    % varargout{1} = photoCurrentTimeAxis;
    return;
else
    varargout{1} = [];
end

%% parse inputs
p = inputParser;
p.addRequired('oi',@isstruct);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('newNoise', true, @islogical);
p.addParameter('emPath', [], @isnumeric);

p.parse(oi,varargin{:});

oi          = p.Results.oi;
currentFlag = p.Results.currentFlag;
newNoise    = p.Results.newNoise;

%% set eye movement path
if isempty(p.Results.emPath)
    assert(~append || isempty(obj.absorptions), ...
        'emPath required when in increment mode');
    emPath = obj.emPositions;
else
    emPath = p.Results.emPath;
    if isempty(obj.absorptions), obj.emPositions = emPath;
    else
        obj.emPositions = [obj.emPositions; emPath];
    end
end

%% extend sensor size
padRows = max(abs(emPath(:, 2)));
padCols = max(abs(emPath(:, 1)));

% We need a copy of the object because ...
cpObj = obj.copy();

% Perhaps because of eye movements?
cpObj.pattern = zeros(obj.rows+2*padRows, obj.cols+2*padCols);

% compute full LMS noise free absorptions
LMS = cpObj.computeSingleFrame(oi, 'fullLMS', true);

% deal with eye movements
absorptions = obj.applyEMPath(LMS, 'emPath', emPath);

% Add photon noise to the whole volume
if obj.noiseFlag
    if (isa(obj, 'coneMosaicHex'))
        % photonNoise is expensive, so only call photonNoise on the
        % non-null cones for a hex mosaic.
        nonNullConeIndices = find(obj.pattern > 1);
        timeSamples = size(absorptions,3);
        absorptions = reshape(permute(absorptions, [3 1 2]), [timeSamples size(obj.pattern,1)*size(obj.pattern,2)]);
        absorptionsCopy = absorptions;
        absorptions = absorptions(:, nonNullConeIndices);
        absorptionsCopy(:, nonNullConeIndices) = obj.photonNoise(absorptions, 'newNoise', newNoise);
        absorptions = permute(reshape(absorptionsCopy, [timeSamples size(obj.pattern,1) size(obj.pattern,2)]), [2 3 1]);
        clear 'absorptionsCopy'
    else
        absorptions = obj.photonNoise(absorptions, 'newNoise', newNoise);
    end
end

% Setting obj.absorptions seems to set the time absorptionsTimeAxis also,
% but incorrectly.  Fix!
% if append
%     obj.absorptions = cat(3, obj.absorptions, absorptions);
% else
%     obj.absorptions = absorptions;
% end

%% If we want the photo current, use the os model

% N.B. You can always calculate the photocurrent later using
%
%   coneMosaic.computeCurrent;
%
% That is BW's preferred method.

current         = [];
% currentTimeAxis = [];  % Used by NC???  Not needed here.
if currentFlag

    if size(obj.absorptions,3) == 1
        disp('Absorptions are a single frame.  No current to calculate.')        
        return;
    else
        % compute the os time axes
        % BW doesn't think this is needed any more
        %
        %         absorptionsTimeAxis = obj.absorptionsTimeAxis;
        %
        %         % Current time axis
        %         dtOS = obj.os.timeStep;
        %         currentTimeAxis = absorptionsTimeAxis(1): dtOS :absorptionsTimeAxis(end);
        
        % Should append be true or false or what?
        obj.current = obj.os.osCompute(cMosaic,'append', append);
        
    end
    
end

end

