
function a = coneAbsorptions(obj,varargin)
% cMosaic.coneAbsorptions('L');
%
% Get absorptions from one the cone classes, possibly in a
% specific ROI
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('coneType','all',@ischar);
p.addParameter('roi',[],@isvector);   % A rect in the array
p.parse(obj,varargin{:});

if ~isempty(p.Results.roi)
    % Clip out the pattern and absorptions
end

% Pull out the relevant cone type data
switch lower(p.Results.coneType)
    case 'l'
        cType = 2;
    case 'm'
        cType = 3;
    case 's'
        cType = 4;
    otherwise
        error('Unknown cone type %s\n',p.Results.coneType);
end

% Should be a better way to do this.  Help, someone.
idx = find(obj.pattern == cType);
[iRows,iCols] = ind2sub(size(obj.pattern), idx);
a = zeros(size(iRows));
for ii=1:length(iRows)
    a(ii) = obj.absorptions(iRows(ii),iCols(ii));
end
end