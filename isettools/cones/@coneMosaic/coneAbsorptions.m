function a = coneAbsorptions(obj, varargin)
% Get isomerizations from one of the cone classes
%
% Sytntax:
%   absorptions = coneAbsorptions(obj, [varargin])
%
% Description:
%    Get the isomerizations from one of the cone classes.
%
% Inputs:
%    obj          - The cone mosaic object
%
% Outputs:
%   'absorptions' - column vector of absorptions in the integration time
%
% Optional key/value pairs:
%    'coneType'   - One of 'L', 'M', or 'S' cone types. Default 'L'
%
% Notes:
%    * [Note: XXX - Maybe we should allow multiple cone types and return a
%      cell array.]
%    * [Note: JNM - Example is not working. Both coneAbsorptions return the
%      error 'index exceeds matrix dimensions'.]

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/22/18  jnm  Formatting

% Examples:
%{
    cm = coneMosaic();
   absorptions = coneAbsorptions(cm)
   absorptions = coneAbsorptions(cm, 'coneType', 'L')
%}

%% Which cone types

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('coneType', 'L', @ischar);
p.parse(obj, varargin{:});

%% Extract the cone type data

switch lower(p.Results.coneType)
    case 'l'
        cType = 2;
    case 'm'
        cType = 3;
    case 's'
        cType = 4;
    otherwise
        error('Unknown cone type %s\n', p.Results.coneType);
end

% Positions of this type
lst = (obj.pattern == cType);

% Absortpions at those positions
a = obj.absorptions(lst);

% May not be necessary
a = a(:);

end