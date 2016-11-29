function [absorptions, current] = compute(obj, oi, varargin)
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
    [absorptions, current] = obj.computeForOISequence(oi,varargin{:});
    return;
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

obj.absorptions = [];
obj.current = [];

%% set eye movement path
if isempty(p.Results.emPath)
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
% vcNewGraphWin; imagesc(absorptions);

% Add photon noise to the whole volume
if obj.noiseFlag
    if (isa(obj, 'coneMosaicHex'))
        % Only call photonNoise on the non-null cones for a hex mosaic.
        nonNullConeIndices = find(obj.pattern > 1);
        timeSamples = size(absorptions,3);
        absorptions = reshape(permute(absorptions, [3 1 2]), [timeSamples size(obj.pattern,1)*size(obj.pattern,2)]);
        absorptionsCopy = absorptions;
        absorptions = absorptions(:, nonNullConeIndices);
        absorptionsCopy(:, nonNullConeIndices) = obj.photonNoise(absorptions, 'newNoise', newNoise);
        absorptions = permute(reshape(absorptionsCopy, [timeSamples size(obj.pattern,1) size(obj.pattern,2)]), [2 3 1]);
        clear 'absorptionsCopy'
    else % Rectangular mosaic
        % Noisy absorptions.  Notice, this does not set the absorptions in
        % the object yet.  We do that below.
        absorptions = obj.photonNoise(absorptions,'newNoise', newNoise);
        % vcNewGraphWin; imagesc(absorptions);
    end
end

% In the single sample case, we set the absorptions in the object.
obj.absorptions = absorptions;

%% If we want the photo current, use the os model

% We recommend that you calculate the photocurrent later using
%   coneMosaic.computeCurrent;
% rather than by setting this flag.

current         = [];
if currentFlag
    if size(obj.absorptions,3) == 1
        disp('Absorptions are a single frame.  No current to calculate.')        
        return;
    else
        % Should append be true or false or what?
        obj.current = obj.os.osCompute(cMosaic);
    end
end

end

