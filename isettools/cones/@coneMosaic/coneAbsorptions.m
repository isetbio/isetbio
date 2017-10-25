function a = coneAbsorptions(obj,varargin)
%CONEABSORPTIONS  Get isomerizations from one of the cone classes
%
%  absorptions = CONEABSORPTIONS(obj,varargin)
%
% Optional parameter name/value pairs chosen from the following:
%   'coneType'  - one of 'L', 'M', or 'S' (default 'L')
%
% Return
%   'absorptions' - column vector of absorptions in the integration
%   time
%
% Example: 
%   absorptions = CONEABSORPTIONS('coneType','L');
%

% HJ ISETBIO Team 2016
%
% Note:  Maybe we should allow multiple cone types and return a cell
% array.
%

%% Which cone types

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('coneType','L',@ischar);
p.parse(obj,varargin{:});

%% Extract the cone type data

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

% Positions of this type
lst = (obj.pattern == cType);

% Absortpions at those positions
a = obj.absorptions(lst);

% May not be necessary
a = a(:);

end