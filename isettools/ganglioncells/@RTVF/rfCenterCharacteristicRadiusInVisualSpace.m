function [anatomicalRFcenterRcDegs, visualRFcenterRcDegs, visualRFcenterConeMap, rotationDegs] = ...
    rfCenterCharacteristicRadiusInVisualSpace(obj, varargin)

    p = inputParser;
    p.addParameter('computeAnatomicalRFcenterRc', true, @islogical);
    p.addParameter('computeVisualRFcenterRc', true, @islogical);
    p.parse(varargin{:});
    computeAnatomicalRFcenterRc = p.Results.computeAnatomicalRFcenterRc;
    computeVisualRFcenterRc = p.Results.computeVisualRFcenterRc;

    theConeIndices = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;
    theConeWeights = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;
    theConeApertureAreas = obj.coneMosaic.computeApertureAreasMetersSquared(theConeIndices) * 1e16;
    theEffectiveOSlengthAttenuationFactors = obj.coneMosaic.computeEffectiveOSlengthAttenuationFactors(theConeIndices);
    
    % Re-center the input cone positions
    theConePositionsDegs = obj.coneMosaic.coneRFpositionsDegs(theConeIndices,:);
    meanConePositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, meanConePositionDegs);

    % Compute the retinal cone input map for the RF center
    retinalRFcenterConeMap = obj.retinalSubregionMapFromPooledConeInputs(...
        theConePositionsDegs, ...
        theConeApertureAreas, ...
        theEffectiveOSlengthAttenuationFactors, ...
        theConeWeights);


    if (~computeAnatomicalRFcenterRc)
        anatomicalRFcenterRcDegs = [];
    else
        % Compute the anatomical RF center Rc
        if (numel(theConeIndices) == 2)
               cone1RFpos = obj.coneMosaic.coneRFpositionsDegs(theConeIndices(1),:);
               cone2RFpos = obj.coneMosaic.coneRFpositionsDegs(theConeIndices(2),:);
               deltaY = cone2RFpos(2)-cone1RFpos(2);
               deltaX = cone2RFpos(1)-cone1RFpos(1);
               forcedOrientationDegs = -atan2d(deltaY, deltaX);
            else
               forcedOrientationDegs = [];
        end
    
        if (numel(theConeIndices) == 1)
            flatTopGaussian = false;
        else
            flatTopGaussian = true;
        end
    
        spatialSupportDegsX = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
        spatialSupportDegsY = obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs;
        theFittedGaussian = RTVF.fitGaussianEllipsoid(...
            spatialSupportDegsX, spatialSupportDegsY, retinalRFcenterConeMap, ...
            'flatTopGaussian', flatTopGaussian, ...
            'forcedOrientationDegs', forcedOrientationDegs, ...
            'globalSearch', false);
        anatomicalRFcenterRcDegs = min(theFittedGaussian.characteristicRadii);
    end


    % Convolve with the LM-weighted PSF
    visualRFcenterConeMap = conv2(retinalRFcenterConeMap, obj.spectrallyWeightedPSFData.LMconeWeighted, 'same');


    % Since RF parameters by Croner&Kaplan were based on optimal
    % orientation gratings (highest resolving SF), we first rotate the
    % visualRFcenterConeMap so that it is narrowest along the x-axis,
    % then sum along the y-dimension, and finally fit a line-weighting function 
    % for a Gaussian 

    % Find best rotation and rotate according to it
    debugRadonTransform = false;
    [visualRFcenterConeMapRotated, rotationDegs] = ...
        RTVF.bestHorizontalResolutionRFmap(visualRFcenterConeMap, [], debugRadonTransform);

    if (~computeVisualRFcenterRc)
        visualRFcenterRcDegs = [];
    else

        % Integrate along y of the 2-D RF
        spatialSampleSize = spatialSupportDegsY(2)-spatialSupportDegsY(1);
        visualRFcenterProfileX = sum(visualRFcenterConeMapRotated,1)*spatialSampleSize;
    
        % Finally, fit a 1D Gaussian line weighting function to the 1D profile 
        spatialSupportXdegs = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
        theFittedGaussianLineWeightingFunction = ...
            RTVF.fitGaussianLineWeightingFunction(spatialSupportXdegs, visualRFcenterProfileX);
        
        % Return the visualRF center Rc (degs)
        visualRFcenterRcDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;

        debugFitting = true;
        if (debugFitting)
            hFig = figure(1); clf;
            set(hFig, 'Position', [10 10 950 950])
            ax = subplot(2,2,1);
            imagesc(...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
                retinalRFcenterConeMap);
            axis 'xy'; axis 'image'
            colormap(gray(1024))
            title('retinal rf center');
        
            ax = subplot(2,2,2);
            imagesc(...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
                visualRFcenterConeMap);
            axis 'xy'; axis 'image'
            title('visual rf center');
    
            ax = subplot(2,2,3);
            imagesc(...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), ...
                obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs + meanConePositionDegs(2), ...
                visualRFcenterConeMapRotated);
            axis 'xy'; axis 'image'
            title(sprintf('visual rf center (rot = %2.1f degs)', rotationDegs));
        
            ax = subplot(2,2,4);
            plot(obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs + meanConePositionDegs(1), visualRFcenterProfileX, 'ks');
            hold on;
            plot(spatialSupportXdegs+ meanConePositionDegs(1), theFittedGaussianLineWeightingFunction.profile, 'r-');
            plot([0 visualRFcenterRcDegs]+ meanConePositionDegs(1), max(theFittedGaussianLineWeightingFunction.profile)*exp(-1)*[1 1], 'b--');
            axis 'square'
        end
    end


end 

