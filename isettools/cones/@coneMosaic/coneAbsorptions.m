function a = coneAbsorptions(obj,varargin)
%CONEABSORPTIONS  Get isomerizations from one the cone classes
%   a = CONEABSORPTIONS(obj,varargin)
%
%   The isomerizations are returned in a column vector.
%
%   Optional parameter name/value pairs chosen from the following:
%
%   'coneType'        'L', 'M', or 'S' (default 'L')
%
%    Example: a = CONEABSORPTIONS('coneType','L');
 
% HJ ISETBIO Team 2016

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('coneType','L',@ischar);
p.parse(obj,varargin{:});

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
a = zeros(size(iRows),1);
for ii=1:length(iRows)
    a(ii) = obj.absorptions(iRows(ii),iCols(ii));
end
end