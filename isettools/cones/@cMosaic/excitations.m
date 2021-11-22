function [roiExcitations, allExcitations] = excitations(cm,varargin)
% Return the excitations within an roi of the cMosaic
%
% Synopsis
%   cMosaic.excitations(roi,varargin);
%   cMosaic.excitations('roi',roi,'all excitations',allE);
%   [roiE, allE] = cMosaic.excitations('roi',roi,'oi',oi);
%
% Brief description
%   Method to get cone excitations from an ROI
%
% Inputs
%   cm - cMosaic object
%
% Optional key/val pairs
%   roi  -  Region of interest
%   allExcitations - Pre-computed excitations
%   oi  - Optical image
%   cone type - {L,M,S}
%
% Outputs
%   roiExcitations - The excitations within the ROI
%   allExcitations - Excitations of the entire mosaic
%
% See also
%   cMosaic.compute(oi);
%

%% Parse
varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cm',@(x)(isa(x,'cMosaic')));
p.addParameter('roi',[],@(x)(isa(x,'regionOfInterest')));
p.addParameter('allexcitations',[],@isnumeric);
p.addParameter('oi',[],@(x)(isstruct(x) && isequal(x.type,'opticalimage')));
p.addParameter('conetype','',@ischar);

p.parse(cm,varargin{:});

roi      = p.Results.roi;
oi       = p.Results.oi;
coneType = p.Results.conetype;

allExcitations = p.Results.allexcitations;
if isempty(allExcitations)
    allExcitations = cm.compute(oi);
end

if isempty(roi)
    roiExcitations = allExcitations;
    return;
end


%% Find indices 

% Here are all the indices of the cones whose positions are within the ROI
idx = roi.indicesOfPointsInside(cm.coneRFpositionsDegs);

% Restrict to a cone type 
switch ieParamFormat(coneType)
    case 'l'
        in = ismember(idx,cm.lConeIndices);
        idx = idx(in);
    case 'm'        
        in = ismember(idx,cm.mConeIndices);
        idx = idx(in);
    case 's'
        in = ismember(idx,cm.sConeIndices);
        idx = idx(in);
    case ''
        % All assumed
    otherwise
        error('Unknown cone type %s\n',coneType);
end

roiExcitations = allExcitations(idx);

end

% in = ismember(idx,cm.lConeIndices);
% idxL = idx(in);
% mean(allExcitations(idxL))

