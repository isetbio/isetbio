function [roiExcitations, roiIdx, allExcitations] = excitations(cm,varargin)
% Return the excitations within a cMosaic ROI
%
%  N.B.  This method will likely be deprecated as we store the
%  excitations within the cMosaic object.  In the future we will be
%  using something like
%
%     cMosaic.get('excitations', varargin)
%
%  as a replacement
%
% Synopsis
%   cMosaic.excitations(roi,varargin);
%   [roiE, roiIdx] = cMosaic.excitations('roi',roi,'all excitations',allE);
%   [roiE, roiIdx, allE] = cMosaic.excitations('roi',roi,'oi',oi);
%
% Brief description
%   Method to read cone excitations from an ROI.  Typically, the
%   excitations are already computed and included as a parameter
%   (allE).  The returns are the excitations within the roi as well as
%   the indices (idx) used to extract data from allE.
%
% Inputs
%   cm - cMosaic object
%
% Optional key/val pairs
%   roi  -  Region of interest
%   allExcitations - Pre-computed excitations
%   oi  - Optical image
%   cone type - {L,M,S}
%   visualize -  Show the ROI overlaid on all the excitations
%
% Outputs
%   roiExcitations - The excitations within the ROI
%   allExcitations - Excitations of the entire mosaic
%   roiIdx         - Indices of the ROI w.r.t allExcitations
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
p.addParameter('visualize',false,@islogical);

p.parse(cm,varargin{:});

roi      = p.Results.roi;
oi       = p.Results.oi;
coneType = p.Results.conetype;

allExcitations = p.Results.allexcitations;
if isempty(allExcitations)
    allExcitations = cm.compute(oi);
end

if isempty(roi)
    % This is no different than cMosaic.compute;
    roiExcitations = allExcitations;
    roiIdx = 1:size(allExcitations,3);
    return;
end


%% Find indices 

% Here are all the indices of the cones whose positions are within the ROI
switch(roi.shape)
    case 'line'
        idx = roi.indicesOfPointsInside(cm.coneRFpositionsDegs);
        %{
        samplingPoints = 500; % sample the perimeter using 1000 points
        pointsPerSample = 10;  % up to 30 points for each sample along the perimeter
        maxDistance = 0.2;     % up to 0.5 units aray from the closest point on the perimeter
        idx = roi.indicesOfPointsAround(cm.coneRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
        % idx = roi.indicesOfPointsAround(cm.coneRFpositionsDegs);
        %}
    otherwise
        idx = roi.indicesOfPointsInside(cm.coneRFpositionsDegs);
end

% Restrict to a cone type 
switch ieParamFormat(coneType)
    case 'l'
        in = ismember(idx,cm.lConeIndices);
        roiIdx = idx(in);
    case 'm'        
        in = ismember(idx,cm.mConeIndices);
        roiIdx = idx(in);
    case 's'
        in = ismember(idx,cm.sConeIndices);
        roiIdx = idx(in);
    case ''
        % All assumed
        roiIdx = idx;
    otherwise
        error('Unknown cone type %s\n',coneType);
end

roiExcitations = allExcitations(roiIdx);

%% Show the ROI if requested

if p.Results.visualize
    cm.plot('roi',allExcitations,'roi',roi);
end

end

