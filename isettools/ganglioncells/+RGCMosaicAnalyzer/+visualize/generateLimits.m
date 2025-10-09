%
% RGCMosaicAnalyzer.visualize.generateLimits()
%

function [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
    domainVisualizationLimits, domainVisualizationTicks, ...
    domainVisualizationLimitsSingleRF, domainVisualizationTicksSingleRF, ...
    sizePixels, sigmaPixels, sizePixelsDegs, sigmaPixelsDegs] = ...
    generateLimits(theMRGCMosaic, theRGCpositionDegs, varargin)

    p = inputParser;
    p.addParameter('spatialSupportTickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.parse(varargin{:});
    spatialSupportTickSeparationArcMin = p.Results.spatialSupportTickSeparationArcMin;

    % The spatial support
    spatialSupportCenterDegs = round(theRGCpositionDegs*10)/10;
    ecc = sqrt(sum(theMRGCMosaic.eccentricityDegs.^2,2));

    if (isempty(spatialSupportTickSeparationArcMin))
        if (ecc > 30)
            spatialSupportTickSeparationArcMin = 60.0;
        elseif (ecc > 20)
            spatialSupportTickSeparationArcMin = 48.0;
        elseif (ecc > 10)
            spatialSupportTickSeparationArcMin = 36.0;
        elseif (ecc > 5)
            spatialSupportTickSeparationArcMin = 24.0;
        elseif (ecc > 3)
            spatialSupportTickSeparationArcMin = 18.0;
        elseif (ecc > 2)
            spatialSupportTickSeparationArcMin = 12.0;
        elseif (ecc > 1)
            spatialSupportTickSeparationArcMin = 9;
        else
            spatialSupportTickSeparationArcMin = 6;
        end
    end

    if (ecc > 30)
        scaleBarMicrons = 40;
    elseif (ecc > 20)
        scaleBarMicrons = 25;
    elseif (ecc > 10)
        scaleBarMicrons = 10;
    elseif (ecc > 5)
        scaleBarMicrons = 5;
    elseif (ecc > 3)
        scaleBarMicrons = 5;
    elseif (ecc > 2)
        scaleBarMicrons = 5;
    elseif (ecc > 1)
        scaleBarMicrons = 5;
    else
        scaleBarMicrons = 2;
    end


    scaleBarDegs = scaleBarMicrons / theMRGCMosaic.inputConeMosaic.micronsPerDegree;

    
    % No scale bar
    % scaleBarDegs = 0;


    domainVisualizationLimits(1:2) = theRGCpositionDegs(1) + [-3 3]*spatialSupportTickSeparationArcMin/60;
    domainVisualizationLimits(3:4) = theRGCpositionDegs(2) + [-3 3]*spatialSupportTickSeparationArcMin/60;

    domainVisualizationTicks = struct(...
        'x', spatialSupportCenterDegs(1) + (-40:40)*spatialSupportTickSeparationArcMin/60, ...
        'y', spatialSupportCenterDegs(2) + (-40:40)*spatialSupportTickSeparationArcMin/60);

    domainVisualizationLimitsSingleRF(1:2) = theRGCpositionDegs(1) + [-3 3]*(1/4)*spatialSupportTickSeparationArcMin/60;
    domainVisualizationLimitsSingleRF(3:4) = theRGCpositionDegs(2) + [-3 3]*(1/4)*spatialSupportTickSeparationArcMin/60;

    domainVisualizationTicksSingleRF = struct(...
        'x', spatialSupportCenterDegs(1) + (-40:40)*(1/4)*spatialSupportTickSeparationArcMin/60, ...
        'y', spatialSupportCenterDegs(2) + (-40:40)*(1/4)*spatialSupportTickSeparationArcMin/60);

end