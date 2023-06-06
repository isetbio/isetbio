function computeInputConeMosaicVisuallyProjectedRcDegs(obj, ...
    responsesFileName, varargin)

    p = inputParser;
    p.addParameter('useParfor', false, @islogical);
    p.addParameter('visualizedResponses', false, @islogical);
    p.addParameter('opticsToEmploy', 'native', @(x)(ismember(x, {'native', 'custom'})));

    p.parse(varargin{:});
    useParfor = p.Results.useParfor;
    visualizeResponses = p.Results.visualizedResponses;
    opticsToEmploy = p.Results.opticsToEmploy;

    % Retrieve the input cone mosaic
    theInputConeMosaic = obj.theRGCMosaic.inputConeMosaic;

    % Load the previously computed inputConeMosaic STF responses
    switch (opticsToEmploy)
        case 'native'
             load(responsesFileName, 'theNativeOpticsParams', ...
                'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
                 'orientationsTested', 'spatialFrequenciesTested', ...
                 'spatialPhasesDegs', 'coneContrasts');

        case 'custom'
             load(responsesFileName, 'theCustomOpticsParams', ...
                'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
                 'orientationsTested', 'spatialFrequenciesTested', ...
                 'spatialPhasesDegs', 'coneContrasts');
        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end

    
    % Transform cone excitation responses to cone modulation responses
    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses== 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);

    % Compute cone mosaic modulation responses
    theConeMosaicModulationSTFresponses = 0 * theConeMosaicSTFresponses;
    for iOri = 1:numel(orientationsTested)
        theConeMosaicModulationSTFresponses(iOri,:,:,:) = ...
            bsxfun(@times, bsxfun(@minus, squeeze(theConeMosaicSTFresponses(iOri,:,:,:)), theConeMosaicNullResponses), ...
                 normalizingResponses);
    end
    clear 'theConeMosaicSTFresponses'

    % Only analyze cones around the horizontal meridian
    theAnalyzedConesROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'line', ...
            'from', [theInputConeMosaic.eccentricityDegs(1)-0.5*theInputConeMosaic.sizeDegs(1) theInputConeMosaic.eccentricityDegs(2)], ...
            'to',   [theInputConeMosaic.eccentricityDegs(1)+0.5*theInputConeMosaic.sizeDegs(1) theInputConeMosaic.eccentricityDegs(2)], ...
            'thickness', 0.05 ...
        ));

    % Find the indices of the cones whose position is within theAnalyzedConesROI 
    analyzedConeIndices = theAnalyzedConesROI.indicesOfPointsInside(theInputConeMosaic.coneRFpositionsDegs);

    % Sort the visualized cone indices according to their x-ecc
    [~,idx] = sort(squeeze(theInputConeMosaic.coneRFpositionsDegs(analyzedConeIndices,1)), 'ascend');
    analyzedConeIndices = analyzedConeIndices(idx);

    % We are visualizing only 32 cones so make sure we cover the entire range
    ii = (1:32)*numel(analyzedConeIndices)/32;
    analyzedConeIndices = analyzedConeIndices(round(ii));

    lmConeIndices = find(...
        (theInputConeMosaic.coneTypes(analyzedConeIndices) == cMosaic.LCONE_ID) | ...
        (theInputConeMosaic.coneTypes(analyzedConeIndices) == cMosaic.MCONE_ID));

    analyzedConeIndices = analyzedConeIndices(lmConeIndices);

    anatomicalRcDegs = zeros(1, numel(analyzedConeIndices));
    visualRcDegs = zeros(1, numel(analyzedConeIndices));
    temporalEquivalentEccDegs = zeros(1,numel(analyzedConeIndices));

    for idx = 1:numel(analyzedConeIndices)
        iCone = analyzedConeIndices(idx);

        % The temporal equivalent ecc for this cone
        tempEquivEccDegs = obj.theRGCMosaic.temporalEquivalentEccentricityForEccentricity(theInputConeMosaic.coneRFpositionsDegs(iCone,:));
        temporalEquivalentEccDegs(idx) = sqrt(sum(tempEquivEccDegs.^2,2));

        % Retrieve the visual STFs for all orientations for this cone
        theConeResponsesAcrossAllOrientationsAndSpatialFrequencies = theConeMosaicModulationSTFresponses(:,:,:,iCone);
    
        [theConeVisualSTF,theSTFsAcrossAllOrientations] = MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
               orientationsTested, spatialFrequenciesTested, ...
               theConeResponsesAcrossAllOrientationsAndSpatialFrequencies);

        fprintf('Computing estimate of the visual Rc for cone %d of %d \n', ...
            idx, numel(analyzedConeIndices));

        % An estimate of the anatomical RcDegs for this cone
        anatomicalRcDegs(idx) = theInputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                           theInputConeMosaic.coneApertureDiametersDegs(iCone);

        rangeForRc = anatomicalRcDegs(idx)*[1 5 20];

        multiStartsNumDoGFit = 256;
        % Fit the visual STF with a DoG model
        [fittedParamsStruct, theFittedSTF] = MosaicPoolingOptimizer.fitGaussianToSubregionSTF(...
                      spatialFrequenciesTested, ...
                      theConeVisualSTF, ...
                      anatomicalRcDegs(idx), ...
                      rangeForRc, ...
                      multiStartsNumDoGFit);

        figure(idx); clf;
        plot(spatialFrequenciesTested, theConeVisualSTF, 'rs', 'LineWidth', 1.5);
        hold on;
        plot(theFittedSTF.sfHiRes, theFittedSTF.subregionSTFHiRes , 'r-', 'LineWidth', 1.5);
        set(gca, 'XScale', 'log')
        
        % The visually projected Rc for this cone
        visualRcDegs(idx) = fittedParamsStruct.finalValues(2);
    end


    figure(99);
    plot(temporalEquivalentEccDegs, visualRcDegs, 'ro-', 'LineWidth', 1.5);
    hold on
    plot(temporalEquivalentEccDegs, anatomicalRcDegs, 'k.-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'XLim', [0.05 30], 'XTick', [0.03 0.1 0.3 1 3 10 30], 'YLim', [0 0.1], 'YTick', 0:0.01:0.1);
    xlabel('temporal equivalent eccentricity (degs)');
    ylabel('Rc (degs)');
    legend({'single cone visual Rc (gaussian STF fit)', 'anatomical Rc'});
end
