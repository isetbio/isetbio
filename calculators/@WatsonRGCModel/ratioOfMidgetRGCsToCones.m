function [ratios, coneRFDensities, mRGCRFDensities] = ratioOfMidgetRGCsToCones(obj, eccXYposDegs, whichEye)
% Compute ratios of midget RGC receptive fields to cones at the requested (x,y) retinal eccentricities (degs)
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   xPos = 0:0.1:1.0;
%   yPos = xPos * 2;
%   eccXYDegs = [xPos(:) yPos(:)];
%   mRGCRFtoConesRatios = WatsonRGCCalc.ratioOfMidgetRGCsToCones(eccXYDegs, 'left');
%
% Description:
%   Method to return the ratio of midget RGC receptive fields to cones at
%   the requested (x,y) retinal eccentricities (specified in degrees,
%   as an [Nx2] matrix).
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    eccXYposDegs              - (x,y) eccentricities ([N x 2] matrix at which to compute the ratios
%    whichEye                  - 'Left' or 'right' eye
%
% Outputs:
%    val                       - Ratio of midget RGC receptive fields to
%                                cones at the requested (x,y) retinal eccentricities 
%                                and eye
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.

% Parse input
    p = inputParser;
    p.addRequired('eccXYposDegs',  ...
        @(x)(isnumeric(x) && (ndims(x) == 2) && (size(x,2) == 2)));
    p.addRequired('whichEye', ...
        @(x)(ischar(x) && ismember(x, {'left', 'right'})));
    p.parse(eccXYposDegs, whichEye);
            
    eccXYposDegs = p.Results.eccXYposDegs;
    whichEye = lower(p.Results.whichEye);
    
    % Compute the on-axis ratios at the requested eccentricity magnitudes
    [onAxisRatios, onAxisAngles, ...
        onAxisMidgetRGCRFDensities, ...
        onAxisConeRFDensities] = computeOnAxisRatios(obj,eccXYposDegs,whichEye);
    
    % Interpolate at the requested eccentricity angles
    ratios = 0.5*interpolateRatios(onAxisRatios,onAxisAngles, eccXYposDegs);
    coneRFDensities = interpolateRatios(onAxisConeRFDensities,onAxisAngles, eccXYposDegs);
    mRGCRFDensities = interpolateRatios(onAxisMidgetRGCRFDensities,onAxisAngles, eccXYposDegs);
end

function ratios = interpolateRatios(onAxisRatios, onAxisAngles, eccXYposDegs)

   eccX = squeeze(eccXYposDegs(:,1));
   eccY = squeeze(eccXYposDegs(:,2));
   requestedAngles = atan2d(eccY,eccX);
   
   % Make sure all angles are > 0
   idx = find(requestedAngles<0);
   requestedAngles(idx) = requestedAngles(idx) + 360;
   
   parfor aa = 1:length(requestedAngles)
        ratios(aa) = interp1(onAxisAngles, onAxisRatios(:,aa), requestedAngles(aa), 'linear');
   end
                    
end

function [onAxisRatios, onAxisAngles, onAxisMidgetRGCRFDensities, onAxisConeRFDensities] = computeOnAxisRatios(obj,eccXYposDegs, whichEye)
    ecc = sqrt(sum(eccXYposDegs.^2,2));

    % Ratios along the 2 horizontal meridians
    meridianName = 'nasal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(ecc', meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(ecc', meridianName, 'Cones per deg2');
    qIndex = 1;
    onAxisRatios(qIndex,:) = midgetRGCRFDensity./coneRFDensity;
    onAxisMidgetRGCRFDensities(qIndex,:) = midgetRGCRFDensity;
    onAxisConeRFDensities(qIndex,:) =  coneRFDensity;
    onAxisAngles(qIndex) = 0;

    meridianName = 'temporal meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(ecc', meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(ecc', meridianName, 'Cones per deg2');
    qIndex = 3;
    onAxisRatios(qIndex,:) = midgetRGCRFDensity./coneRFDensity;
    onAxisMidgetRGCRFDensities(qIndex,:) = midgetRGCRFDensity;
    onAxisConeRFDensities(qIndex,:) =  coneRFDensity;
    onAxisAngles(qIndex) = 180;        
    
    if (strcmp(whichEye, 'right'))
        tmp = onAxisRatios(1,:);
        onAxisRatios(1,:) = onAxisRatios(3,:);
        onAxisRatios(3,:) = tmp;
        tmp = onAxisMidgetRGCRFDensities(1,:);
        onAxisMidgetRGCRFDensities(1,:) = onAxisMidgetRGCRFDensities(3,:);
        onAxisMidgetRGCRFDensities(3,:) = tmp;
        tmp = onAxisConeRFDensities(1,:);
        onAxisConeRFDensities(1,:) = onAxisConeRFDensities(3,:);
        onAxisConeRFDensities(3,:) = tmp;
        clear 'tmp';
    end
    
    % Ratios along the 2 vertical meridians
    meridianName = 'superior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(ecc', meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(ecc', meridianName, 'Cones per deg2');
    qIndex = 2;
    onAxisRatios(qIndex,:) = midgetRGCRFDensity./coneRFDensity;
    onAxisMidgetRGCRFDensities(qIndex,:) = midgetRGCRFDensity;
    onAxisConeRFDensities(qIndex,:) =  coneRFDensity;
    onAxisAngles(qIndex) = 90;
    
    meridianName = 'inferior meridian';
    midgetRGCRFDensity = obj.midgetRGCRFDensity(ecc', meridianName, 'RFs per deg2');
    [~, coneRFDensity] = obj.coneRFSpacingAndDensity(ecc', meridianName, 'Cones per deg2');
    qIndex = 4;
    onAxisRatios(qIndex,:) = midgetRGCRFDensity./coneRFDensity;
    onAxisMidgetRGCRFDensities(qIndex,:) = midgetRGCRFDensity;
    onAxisConeRFDensities(qIndex,:) =  coneRFDensity;
    onAxisAngles(qIndex) = 270;
    
    % Wrap around
    onAxisRatios(5,:) = onAxisRatios(1,:);
    onAxisMidgetRGCRFDensities(5,:) = onAxisMidgetRGCRFDensities(1,:);
    onAxisConeRFDensities(5,:) = onAxisConeRFDensities(1,:);
    onAxisAngles(5) = onAxisAngles(1)+360;
end

    