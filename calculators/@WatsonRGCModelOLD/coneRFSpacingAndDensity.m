function [coneRFSpacing, coneRFDensity, meridianAngle] = coneRFSpacingAndDensity(obj, ecc, meridianName, meridianSpace, whichEye, eccUnits, returnUnits)
% Return cone receptive field spacing and density at the requested meridian
% and eccentricities for the specified eye, meridian space.
%
% Syntax:
%   WatsonRGCCalc = WatsonRGCModel();
%   ecc = 0:0.1:10;
%   meridianName = 'superior meridian';
%   meridianSpace = 'visual';
%   eye = 'left';
%   eccUnits = 'deg';
%   returnUnits = 'Cones per deg2';
%   [coneSpacingMM, coneDensityPerMM2] = WatsonRGCCalc.coneRFSpacingAndDensity(ecc, meridianName, meridianSpace, eye, eccUnits, 'Cones per mm2')
%   [coneSpacingDeg, coneDensityPerDeg2] = WatsonRGCCalc.coneRFSpacingAndDensity(ecc, meridianName, meridianSpace, eye, eccUnits, 'Cones per deg2')
%
% Description:
%   Method to return cone spacing as a function of the requested meridian and 
%   eccentricities, specified in visual degrees (with returned spacing/density specified either in deg or mm).
%
% Inputs:
%    obj                       - The WatsonRGCModel object
%    ecc                       - Eccentricities at which to compute RF densities
%    meridianName              - Meridian for which to compute RF densities
%    meridianSpace             - Space in which meridian is specified ('visual' or 'retinal')
%    whichEye                  - 'left' or 'right'
%    eccUnits                  - Units at which the passed eccentricities are specified
%    returnUnits               - Return units, either 'Cones per mm2'
%                                or 'Cones per deg2'
% Outputs:
%    val                       - Cone spacing at the requested eccentricities
% 
% References:
%    Watson (2014). 'A formula for human RGC receptive field density as
%    a function of visual field location', JOV (2014), 14(7), 1-17.
%
% History:
%    11/11/19  NPC, ISETBIO Team     Wrote it.
    
  
    switch (eccUnits)
        case 'deg'
            % Convert ecc from degs to retinal MMs
            eccMM = obj.rhoDegsToMMs(ecc);
            eccDegs = ecc;
        case 'mm'
            eccMM = ecc;
            eccDegs = obj.rhoMMsToDegs(eccMM);
        otherwise
            error('eccUnits must be either ''deg'' or ''mm''. ''%s'' not recognized', eccUnits);
    end
    
    
    % Call the isetbio function coneSizeReadData to read-in the Curcio '1990
    % cone spacing/density data.
    meridianAngle = obj.meridianAngle(meridianName, meridianSpace, whichEye);
    [~, ~, densityConesPerMM2] = coneSizeReadData('eccentricity', eccMM, ...
                                        'angle', meridianAngle*ones(1,numel(eccMM)), ...
                                        'eccentricityUnits', 'mm', ...
                                        'whichEye', whichEye, ...
                                        'useParfor', false);
        
    % Apply correction for the fact that the isetbio max cone density (18,800 cones/deg^2) 
    % does not agree with Watson's (obj.dc0 =  14,804.6 cones/deg^2), and the fact that if we do not
    % apply this correction we get less than 2 mRGCs/cone at foveal eccentricities. We
    % apply this correction only for ecc <= 0.18 degs
    correctForFovealEcc = true;
    if (correctForFovealEcc)                                
        eccLimit = 0.18;
        
        WatsonModelMaxConeDensityPerDeg2 = obj.dc0;
        [~,~,ISETBioMaxConeDensityPerMM2] = coneSizeReadData('eccentricity', 0, 'angle', 0);
        ISETBioMaxConeDensityPerDeg2 = ISETBioMaxConeDensityPerMM2 * obj.alpha(0);
    
        correctionFactorMax = ISETBioMaxConeDensityPerDeg2 - WatsonModelMaxConeDensityPerDeg2;
        correctionFactorMax = correctionFactorMax / obj.alpha(0);
        
        idx = find(abs(eccDegs)<=eccLimit);
        if (~isempty(idx))
            indicesToBeCorrected = idx;
            correctionFactors = correctionFactorMax.*(eccLimit-eccDegs(indicesToBeCorrected))/eccLimit;
            densityConesPerMM2(indicesToBeCorrected) = densityConesPerMM2(indicesToBeCorrected) - correctionFactors;
        end
    end
    
    % In ConeSizeReadData, spacing is computed as sqrt(1/density). This is
    % true for a rectangular mosaic. For a hex mosaic, spacing = sqrt(2.0/(3*density)).
    spacingMM = sqrt(2.0./(sqrt(3.0)*densityConesPerMM2));
    spacingMeters = spacingMM * 1e-3;
     
    switch (returnUnits)
        case 'cones per mm2'
            coneRFSpacing = spacingMeters * 1e3;
            coneRFDensity = densityConesPerMM2;
            
        case 'cones per deg2'
            spacingMM = spacingMeters * 1e3;
            % Convert cone spacing in mm to cone spacing in degs at all eccentricities
            coneRFSpacing = obj.rhoMMsToDegs(spacingMM+eccMM)-obj.rhoMMsToDegs(eccMM); 
            
            % Convert cone density from per mm2 to per deg2
            % Compute mmSquaredPerDegSquared conversion factor for the
            % eccentricities (ecc specified in degs)
            mmSquaredPerDegSquared = obj.alpha(eccDegs);
            coneRFDensity = densityConesPerMM2 .* mmSquaredPerDegSquared;
        otherwise
            error('Density units must be either ''cones per mm2'' or ''cones per deg2''.');
    end
    
end