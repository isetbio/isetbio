function [absorptions, current, currentTimeAxis, varargout] = compute(obj, oi, varargin)
% Compute the pattern of cone absorptions and typically the
% photocurrent
%    [absorptions, current, currentTimeAxis] = cMosaic.compute(oi);
%
% Inputs:
%   oi  - optical image, or oiSequence.  See oiCreate for more details
%
% Optional inputs:
%   currentFlag  - logical, whether to compute photocurrent
%   newNoise     - logical, whether to use new random seed
%   append       - logical, whether to append to existing data
%   emPath       - eye movement path in nx2 matrix. This
%                  parameter shadows obj.emPositions and is
%                  required when append is true
%
% Outputs:
%   absorptions  - cone photon absorptions
%   current      - cone photocurrent
%
% Notes:
%   If you have absorptions and want to compute photocurrent
%   only, use
%     pRate = absorptions / obj.integrationTime;
%     current = obj.os.osCompute(pRate, obj.pattern);
%
%   When append is true, the stored data will increment.
%   However, the returned absorptions and current are for the
%   current oi only.
%
% HJ ISETBIO Team 2016

%% Check if an oi sequence. 
% Send to the specialized compute in that case.
% Otherwise, just carry on.
if isequal(class(oi),'oiSequence')
    [absorptions, current, currentTimeAxis, photoCurrentTimeAxis] = ...
        obj.computeForOISequence(oi,varargin{:});
    varargout{1} = photoCurrentTimeAxis;
    return;
else
    varargout{1} = [];
end

%% parse inputs
p = inputParser;
p.addRequired('oi',@isstruct);
p.addParameter('currentFlag', false, @islogical);
p.addParameter('newNoise', true, @islogical);
p.addParameter('append', false, @islogical);
p.addParameter('emPath', [], @isnumeric);

p.parse(oi,varargin{:});
oi = p.Results.oi;
currentFlag = p.Results.currentFlag;
newNoise = p.Results.newNoise;
append = p.Results.append;

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
    absorptions = obj.photonNoise(absorptions, ...
        'newNoise', newNoise);
end

% Setting obj.absorptions seems to set the time absorptionsTimeAxis also,
% but incorrectly.  Fix!
if append
    obj.absorptions = cat(3, obj.absorptions, absorptions);
else
    obj.absorptions = absorptions;
end

%% If we want the photo current, use the os model
%
% You can always calculate the photocurrent later using
% coneMosaic..computeCurrent;

current = [];
currentTimeAxis = [];
if currentFlag

    if size(obj.absorptions,3) == 1
        disp('Absorptions are a single frame.  No current to calculate.')        
        return;
    else
        % compute the os time axes
        absorptionsTimeAxis = obj.absorptionsTimeAxis;
        
        % Current time axis
        dtOS = obj.os.timeStep;
        currentTimeAxis = absorptionsTimeAxis(1): dtOS :absorptionsTimeAxis(end);
        
        % Resample the absorptions over time
        %resampledAbsorptionsSequence = ...
        %    coneMosaic.tResample(absorptions, obj.pattern, absorptionsTimeAxis, currentTimeAxis);
        
        % Keep the total number of photon absorptions constant by
        % correcting for the new delta time (dtOS)
        % pRate = resampledAbsorptionsSequence/dtOS;
        
        % Should append be true or false or what?
        current = obj.os.osCompute(cMosaic,'append', append);
        
    end
    
end

end

