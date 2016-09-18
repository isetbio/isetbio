function absorptions = applyEMPath(obj, LMS, varargin)
% apply eye movement path and pick out corresponding cone type
%    absorptions = applyEMPath(obj, LMS, emPath, varargin)
%
% Inputs:
%   obj     - cone mosaic object
%   LMS     - full LMS noise free absorptions
%
% Optional inputs (key-value pairs)
%   padRows - rows padded
%   padCols - columns padded
%   emPath  - eye movement path
%
% Outputs:
%   absorptions - cone absorptions with eye movements
%
% 
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