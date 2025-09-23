%
% RGCMosaicAnalyzer.visualize.generateLimits()
%

function [scaleBarDegs, scaleBarMicrons, spatialSupportTickSeparationArcMin, spatialSupportCenterDegs, ...
    domainVisualizationLimits, domainVisualizationTicks, ...
    domainVisualizationLimitsSingleRF, domainVisualizationTicksSingleRF, ...
    sizePixels, sigmaPixels, sizePixelsDegs, sigmaPixelsDegs] = ...
    generateLimits(theMRGCMosaic, theRGCpositionDegs)

    % The spatial support
    spatialSupportCenterDegs = round(theMRGCMosaic.eccentricityDegs*10)/10;
    ecc = sqrt(sum(spatialSupportCenterDegs.^2,2));

    if (ecc > 30)
        scaleBarDegs = 0.25;
        scaleBarMicrons = 40;
        spatialSupportTickSeparationArcMin = 60.0;
    elseif (ecc > 20)
        scaleBarDegs = 0.1;
        scaleBarMicrons = 25;
        spatialSupportTickSeparationArcMin = 40.0;
    elseif (ecc > 10)
        scaleBarDegs = 0.08;
        scaleBarMicrons = 10;
        spatialSupportTickSeparationArcMin = 20.0;
    elseif (ecc > 5)
        scaleBarDegs = 0.07;
        scaleBarMicrons = 5;
        spatialSupportTickSeparationArcMin = 13.5;
    elseif (ecc > 3)
         scaleBarDegs = 0.07;
        scaleBarMicrons = 5;
        spatialSupportTickSeparationArcMin = 13.5;
    elseif (ecc > 2)
         scaleBarDegs = 0.07;
        scaleBarMicrons = 5;
        spatialSupportTickSeparationArcMin = 13.5;
    elseif (ecc > 1)
         scaleBarDegs = 0.07;
        scaleBarMicrons = 5;
        spatialSupportTickSeparationArcMin = 12;
    else
        scaleBarMicrons = 2;
        spatialSupportTickSeparationArcMin = 6;
        scaleBarDegs = spatialSupportTickSeparationArcMin/60;
    end

    scaleBarDegs = scaleBarMicrons / theMRGCMosaic.inputConeMosaic.micronsPerDegree;

    % No scale bar
    %scaleBarDegs = 0;

    % XY lims
    spatialSupportRangeArcMin = 3 * spatialSupportTickSeparationArcMin;
    maxXY = round(spatialSupportRangeArcMin/2);
    spatialSupportDegs = (-maxXY:0.05:maxXY)/60;
    spatialSupportXYDegs(:,1) = spatialSupportCenterDegs(1) + spatialSupportDegs;
    spatialSupportXYDegs(:,2) = spatialSupportCenterDegs(2) + spatialSupportDegs;
    dx = (spatialSupportDegs(end)-spatialSupportDegs(1))*0.05;
    XLims = spatialSupportCenterDegs(1) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];
    YLims = spatialSupportCenterDegs(2) + [spatialSupportDegs(1)-dx spatialSupportDegs(end)+dx];

    domainVisualizationLimits(1:2) = theRGCpositionDegs(1) + [-3 3]*spatialSupportTickSeparationArcMin/60;
    domainVisualizationLimits(3:4) = theRGCpositionDegs(2) + [-3 3]*spatialSupportTickSeparationArcMin/60;

    domainVisualizationTicks = struct(...
        'x', spatialSupportCenterDegs(1) + (-10:10)*spatialSupportTickSeparationArcMin/60, ...
        'y', spatialSupportCenterDegs(2) + (-10:10)*spatialSupportTickSeparationArcMin/60);

    domainVisualizationLimitsSingleRF(1:2) = theRGCpositionDegs(1) + [-3 3]*(1/4)*spatialSupportTickSeparationArcMin/60;
    domainVisualizationLimitsSingleRF(3:4) = theRGCpositionDegs(2) + [-3 3]*(1/4)*spatialSupportTickSeparationArcMin/60;

    domainVisualizationTicksSingleRF = struct(...
        'x', spatialSupportCenterDegs(1) + (-10:10)*(1/4)*spatialSupportTickSeparationArcMin/60, ...
        'y', spatialSupportCenterDegs(2) + (-10:10)*(1/4)*spatialSupportTickSeparationArcMin/60);
end