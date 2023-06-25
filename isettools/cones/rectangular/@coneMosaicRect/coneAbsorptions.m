function [a, locs, lst] = coneAbsorptions(obj, varargin)
% Get isomerizations from one of the cone classes
%
% Sytntax:
%   [absorptions, locs, lst] = coneAbsorptions(obj, [varargin])
%
% Description:
%    Get the isomerizations from one of the cone classes.
%
% Inputs:
%    obj          - The cone mosaic object
%
% Outputs:
%   'absorptions' - Selected absorptions or all absorptions (vector) 
%   'locs'        - Locations of selected cones (n x 2)
%   'lst'         - Logical indices into coneMosaic.absorptions where the
%                   values were found (row x col)
%
% Optional key/value pairs:
%    'cone type'   - One of 'L', 'M', or 'S' cone types or LMS for all.
%                    Default is LMS
%    'units'       - Spatial units (default is meters)
%
% Notes:
%    * [Note: XXX - Maybe we should allow multiple cone types and return a
%      cell array.]
%
% See also:
% 

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/22/18  jnm  Formatting
%    04/07/18  dhb  Fix example.

% Examples:
%{
   scene = sceneCreate;
   oi = oiCreate('human');
   oi = oiCompute(oi,scene);
   cm = coneMosaic();
   cm.compute(oi);
   absorptions = cm.coneAbsorptions;
   [absorptions, locsL] = coneAbsorptions(cm, 'cone type', 'l', 'units','um');
   [absorptions, locsM] = coneAbsorptions(cm, 'cone type', 'm', 'units','um');
   ieNewGraphWin; 
   plot(locsL(:,1),locsL(:,2),'ro'); hold on;
   plot(locsM(:,1),locsM(:,2),'gx'); 
   axis equal; grid on; xlabel('Position (um)'); ylabel('Position (um)');
%}

%% Which cone types

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('conetype', 'lms', @ischar);
p.addParameter('units', 'm', @ischar);

p.parse(obj, varargin{:});
units = p.Results.units;

% Remains empty for lms case
lst = [];

%% Extract the cone type data

switch lower(p.Results.conetype)
    case 'l'
        cType = 2;
    case 'm'
        cType = 3;
    case 's'
        cType = 4;
    case 'lms'
        a = obj.absorptions;
        a = a(:);
        return;
    otherwise
        error('Unknown cone type %s\n', p.Results.conetype);
end

% Positions of this type
lst = (obj.pattern == cType);

% Absortpions at those positions as a vector
a = obj.absorptions(lst);

locs = obj.coneLocs(lst,:)*ieUnitScaleFactor(units);

end