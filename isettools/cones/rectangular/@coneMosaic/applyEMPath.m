function absorptions = applyEMPath(obj, LMS, varargin)
% Apply eye movement path and return the absorptions.
%
% Syntax:
%   absorptions = applyEMPath(obj, LMS, emPath, [varargin])
%
% Description:
%    Apply the eye movement path and return the absorptions.
%
% Inputs:
%    obj         - A cone mosaic object
%    LMS         - The full LMS noise free absorptions
%
% Outputs:
%    absorptions - The cone absorptions with eye movements
%
% Optional key/value pairs:
%    'name'      - Mosaic name
%    'padRows'   - Rows padded
%    'padCols'   - Columns padded
%    'emPath'    - Eye movement path
%
% Notes:
%    * [Note: DHB - SOME COMMENTS ABOUT HOW THIS EM STUFF WORKS, SOMEWHERE,
%      WITH A POINTER HERE WOULD BE VERY HELPFUL. WHAT ARE THE UNITS OF EM
%      POSITIONS, WHO IS IN CHARGE OF MAKING SURE THE MOSAIC IS BIG ENOUGH,
%      ETC.?]
%
%   See Also:
%    coneMosaic
%

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/22/18  jnm  Formatting

% parse inputs
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addRequired('LMS', @isnumeric);
p.addParameter('padRows', [], @isnumeric);
p.addParameter('padCols', [], @isnumeric);
p.addParameter('emPath', obj.emPositions, @isnumeric);

p.parse(obj, LMS, varargin{:});
emPath = p.Results.emPath;
padRows = p.Results.padRows;
padCols = p.Results.padCols;

xpos = emPath(:, 1);
ypos = emPath(:, 2);
if isempty(padRows), padRows = max(abs(ypos)); end
if isempty(padCols), padCols = max(abs(xpos)); end

% prepare parameters for cone type mask
mask = zeros(obj.rows, obj.cols, 3); % locations for cones
% L, M, S cones
for ii = 2 : 4, mask(:,:,ii - 1) = double(obj.pattern == ii); end

% select out the subframe given the eye position and cone type
absorptions = zeros(obj.rows, obj.cols, size(emPath, 1));
for ii = 1 : size(emPath, 1)
    % sample by position
    cropLMS = LMS((1 + padRows + ypos(ii)):(end - padRows + ypos(ii)), ...
        (1 + padCols - xpos(ii)):(end - padCols - xpos(ii)), :);

    % sample by conetype
    absorptions(:, :, ii) = sum(cropLMS .* mask, 3);
end
end