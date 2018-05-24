function [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = computeFixationMap(...
        timeAxis, emPaths, emPosRange, emPosDelta, varargin)
% Compute a density map of the area spanned by an eye movement path
%
% Syntax:
%   [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
%       fixationMapXSlice, fixationMapYSlice] = computeFixationMap(...
%       timeAxis, emPaths, emPosRange, emPosDelta, [varargin])
%
% Description:
%    The purpose of this function is to compute the fixation map for
%    a fixationalEM object. This will not work if fixationalEM.compute
%    has not already been executed.
%
%    The code below contains examples. To access, type 'edit
%    computeFixationMap.m' into your Command Window.
%
% Inputs:
%    timeAxis             - Vector. Time axis for fixationalEM object.
%    emPaths              - Matrix. The eye movement paths. Format is
%                           Trials x samples x 2.
%    emPosRange           - Vector. The eye movement position ranges.
%                           This is a 1x2 Vector containing min & max.
%    emPosDelta           - Numeric. The eye movement position delta.
%    varargin             - (Optional) Additional arguments, as needed.
%                           See the optional key/values section.
%
% Outputs:
%    fixationMap          - Matrix. The fixation map. It contains
%                           emPosRange * emPosDelta entries per side.
%    fixationMapSupportX  - Vector. The X support for the fixation map.
%                           It is a 1 x emPosRange * emPosDelta vector.
%    fixationMapSupportY  - Vector. The Y support for the fixation map.
%                           It is a 1 x emPosRange * emPosDelta vector.
%    fixationMapXSlice    - Vector. A X slice for the fixation map. It
%                           is a 1 x emPosRange * emPosDelta vector.
%    fixationMapYSlice    - Vector. A Y slice for the fixation map. It
%                           is a emPosRange * emPosDelta x 1 vector.
%
% Optional key/value pairs:
%    'maxDurationSeconds' - Numeric. The duration for which the fixation
%                           map is computed, in seconds.
%                           Default is infinite, which will use the entire
%                           emPath.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments

% Examples:
%{
    fixEMobj = fixationalEM();
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'none';
    fixEMobj.randomSeed = 1;
    fixEMobj.compute(1.0, 1/1000, 2048, false, 'useParfor', true);

    [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = ...
        fixEMobj.computeFixationMap(fixEMobj.timeAxis, ...
        fixEMobj.emPosArcMin, [-20 20], 0.5);
%}

% Parse inputs
p = inputParser;
p.addRequired('timeAxis', @isnumeric);
p.addRequired('emPaths', @isnumeric);
p.addRequired('emPosRange', @isnumeric);
p.addRequired('emPosDelta', @isnumeric);
addParameter(p, 'maxDurationSeconds', Inf, @isnumeric);
p.parse(timeAxis, emPaths, emPosRange, emPosDelta, varargin{:});

% Only analyze span within the maxDurationSeconds
idx = find(timeAxis <= p.Results.maxDurationSeconds);
emPaths = emPaths(:, idx, :);

xEdges = emPosRange(1):emPosDelta:emPosRange(end) + emPosDelta;
yEdges = emPosRange(1):emPosDelta:emPosRange(end) + emPosDelta;
xPos = squeeze(emPaths(:, :, 1));
yPos = squeeze(emPaths(:, :, 2));
fixationMap = histcounts2(xPos(:), yPos(:), xEdges, yEdges);
fixationMap = fixationMap / max(fixationMap(:));
fixationMapSupportX = xEdges(1:end - 1) + emPosDelta / 2;
fixationMapSupportY = yEdges(1:end - 1) + emPosDelta / 2;

[~, midRow] = min(abs(yEdges));
fixationMapXSlice = squeeze(fixationMap(midRow, :));
fixationMapYSlice = squeeze(fixationMap(:, midRow));

end