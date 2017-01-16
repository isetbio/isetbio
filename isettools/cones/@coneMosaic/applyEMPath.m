function absorptions = applyEMPath(obj, LMS, varargin)
%APPLYEMPATH  Apply eye movement path and return absorptoins absorptions.
%    absorptions = APPLYEMPATH(obj, LMS, emPath, varargin)
%
%    Inputs:
%    obj     - cone mosaic object
%    LMS     - full LMS noise free absorptions
%
%    Optional parameter name/value pairs chosen from the following:
%
%    'name'            Mosaic name
%    'padRows'         Rows padded
%    'padCols'         Columns padded
%    'emPath'          Eye movement path
%
%    Outputs:
%    absorptions - cone absorptions with eye movements
%
%    See also CONEMOSAIC.

% HJ ISETBIO Team 2016

% parse inputs
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addRequired('LMS', @isnumeric);
p.addParameter('padRows', [], @isnumeric);
p.addParameter('padCols', [], @isnumeric);
p.addParameter('emPath', obj.emPositions, @isnumeric);

p.parse(obj, LMS, varargin{:});
emPath = p.Results.emPath;
padRows = p.Results.padRows; padCols = p.Results.padCols;

xpos = emPath(:, 1); ypos = emPath(:, 2);
if isempty(padRows), padRows = max(abs(ypos)); end
if isempty(padCols), padCols = max(abs(xpos)); end

% prepare parameters for cone type mask
mask = zeros(obj.rows, obj.cols, 3); % locations for cones
for ii = 2 : 4 % L, M, S
    mask(:,:,ii-1) = double(obj.pattern == ii);
end

% select out the subframe given the eye position and cone type
absorptions = zeros(obj.rows, obj.cols, size(emPath, 1));
for ii = 1 : size(emPath, 1)
    % sample by position
    cropLMS = LMS((1+padRows+ypos(ii)):(end-padRows+ypos(ii)), ...
        (1+padCols-xpos(ii)):(end-padCols-xpos(ii)), :);
    
    % sample by conetype
    absorptions(:, :, ii) = sum(cropLMS .* mask, 3);
end
end