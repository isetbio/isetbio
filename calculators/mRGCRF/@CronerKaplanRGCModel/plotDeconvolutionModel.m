function plotDeconvolutionModel(deconvolutionModel)
    
    metaData = visualizeCenterDeconvolutionModel(deconvolutionModel.center);
    visualizeSurroundDeconvolutionModel(deconvolutionModel.surround, metaData);
    
end

function visualizeSurroundDeconvolutionModel(deconvolutionModel, metaData)
    hFig = figure(2);
    set(hFig, 'Name', sprintf('RF surround deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    assert(metaData.subjectID == deconvolutionModel.subjectID, 'subjectID for center and surround do not match');
    assert(strcmp(metaData.quadrantName, deconvolutionModel.quadrant), 'quadrant names for center and surround do not match');
    
    eccsNum = numel(deconvolutionModel.tabulatedEccentricityRadii);
    conesNum = size(deconvolutionModel.nominalSurroundRetinalCharacteristicRadii,2);
    surroundsNum = size(deconvolutionModel.nominalSurroundRetinalCharacteristicRadii,3);
    nominalRetinalCharacteristicRadii = reshape(deconvolutionModel.nominalSurroundRetinalCharacteristicRadii, ...
        [eccsNum conesNum*surroundsNum]);
    visualCharacteristicRadii = reshape(deconvolutionModel.characteristicRadiusDegs, ...
        [eccsNum conesNum*surroundsNum]);
    peakSensitivities = reshape(deconvolutionModel.peakSensitivity, ...
        [eccsNum conesNum*surroundsNum]);
    

    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.07, ...
            'bottomMargin', 0.1, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.15, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.1, ...
            'rowsNum', 1, ...
            'colsNum', 2);

    
    ax = theAxesGrid{1,1};
    plot(ax, deconvolutionModel.tabulatedEccentricityRadii, nominalRetinalCharacteristicRadii ./ visualCharacteristicRadii, 'ko');
    set(ax, 'YLim', [0 1]);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'retinal:visual radius');
    
    ax = theAxesGrid{1,2};
    plot(ax, deconvolutionModel.tabulatedEccentricityRadii, peakSensitivities, 'ko');
    set(ax, 'YLim', [0 1]);
    xlabel(ax, 'eccentricity (degs)');
    ylabel(ax, 'peak sensitivity');
    pause
    
end

function metaData = visualizeCenterDeconvolutionModel(deconvolutionModel)

    w = WatsonRGCModel();
    retinalNasalEccentricitiesDegs = deconvolutionModel.tabulatedEccentricityRadii;
    coneRFSpacingsDegs = w.coneRFSpacingAndDensityAlongMeridian(retinalNasalEccentricitiesDegs, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureDiameterDegs = WatsonRGCModel.coneApertureToDiameterRatio * coneRFSpacingsDegs;
    coneRadiiDegs = 0.5*coneApertureDiameterDegs;
    

    maxConeInputsNum = 0;
    
    % Ensure all ecc data are for the same subject, quadrant, and sensitivity range
    for eccIndex = 1:numel(deconvolutionModel.metaData)
        metaData = deconvolutionModel.metaData{eccIndex};
        if (eccIndex == 1)
            theMetaData = metaData;
        end
        assert(strcmp(deconvolutionModel.quadrant, metaData.quadrantName), ...
            sprintf('quadrant name mismatch at ecc %2.1f degs: ''%s'' vs ''%s''.', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), theMetaData.quadrantName, metaData.quadrantName));
        
        assert(deconvolutionModel.subjectID == metaData.subjectID, ...
            sprintf('subject ID mismatch at ecc %2.1f degs: ''%s'' vs ''%s''.', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), theMetaData.subjectID, metaData.subjectID));
        
        assert(theMetaData.sensitivityRangeOverWhichToMatchSFtuning(1) == metaData.sensitivityRangeOverWhichToMatchSFtuning(1), ...
            sprintf('sensitivityRangeOverWhichToMatchSFtuning mismatch(1) at ecc %2.1f degs: ''%s'' vs ''%s''.', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), ...
            theMetaData.sensitivityRangeOverWhichToMatchSFtuning(1), ...
            metaData.sensitivityRangeOverWhichToMatchSFtuning(1)));
        
        assert(theMetaData.sensitivityRangeOverWhichToMatchSFtuning(2) == metaData.sensitivityRangeOverWhichToMatchSFtuning(2), ...
            sprintf('sensitivityRangeOverWhichToMatchSFtuning mismatch(1)  at ecc %2.1f degs: ''%s'' vs ''%s''.', ...
            deconvolutionModel.tabulatedEccentricityRadii(eccIndex), ...
            theMetaData.sensitivityRangeOverWhichToMatchSFtuning(2), ...
            metaData.sensitivityRangeOverWhichToMatchSFtuning(2)));   
        
        numberOfConeInputsForEccentricity = sum(~isnan(squeeze(deconvolutionModel.centerConeInputsNum(eccIndex,:))));
        fprintf('Number of cone inputs in RF center examined for ecc %2.1f degs: %d\n', deconvolutionModel.tabulatedEccentricityRadii(eccIndex), ...
            numberOfConeInputsForEccentricity);
        
        if (numberOfConeInputsForEccentricity > maxConeInputsNum)
            maxConeInputsNum = numberOfConeInputsForEccentricity;
        end
        
    end

    ck = CronerKaplanRGCModel('generateAllFigures', false);
    ck.setupPlotLab(0, 25, 15);
    
    hFig = figure(1); clf;
    set(hFig, 'Name', sprintf('RF center deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    lineColors = brewermap(size(deconvolutionModel.centerConeInputsNum,2), '*blues');
    
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.07, ...
            'bottomMargin', 0.1, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.15, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.1, ...
            'rowsNum', 1, ...
            'colsNum', 2);

    
    ax = theAxesGrid{1,1};
    legends = {};
    for inputConesNum = 1:size(deconvolutionModel.centerConeInputsNum,2)
        switch (inputConesNum)
            case 1
                markerType = 'o';
            case 2
                markerType = '+';
            case 3
                markerType = '^';
            case 4
                markerType = 'd';
            case 5
                markerType = 'x';
            case 6
                markerType = '*';
            otherwise
                markerType = '.';
        end
        
        plot(ax,deconvolutionModel.tabulatedEccentricityRadii, deconvolutionModel.characteristicRadiusDegs(:,inputConesNum), ...
            'o-', 'Color', squeeze(lineColors(inputConesNum,:)), 'Marker', markerType);
        hold(ax, 'on');
        legends{inputConesNum} = sprintf('%d-cones', inputConesNum);
    end
    plot(ax,deconvolutionModel.tabulatedEccentricityRadii, coneRadiiDegs/3, 'k--');
    
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'characteristic radius (degs)');

    ax = theAxesGrid{1,2};
    for inputConesNum = 1:size(deconvolutionModel.centerConeInputsNum,2)
        switch (inputConesNum)
            case 1
                markerType = 'o';
            case 2
                markerType = '+';
            case 3
                markerType = '^';
            case 4
                markerType = 'd';
            case 5
                markerType = 'x';
            case 6
                markerType = '*';
            otherwise
                markerType = '.';
        end
        plot(ax,deconvolutionModel.tabulatedEccentricityRadii, deconvolutionModel.peakSensitivity(:,inputConesNum), ...
            'o-', 'Color', squeeze(lineColors(inputConesNum,:)), 'Marker', markerType);
        hold(ax, 'on');
    end
    xlabel(ax,'eccentricity (degs)');
    ylabel(ax,'peak sensitivity');
    set(ax, 'YLim', [0 1]);
    legend(ax, legends, 'Location', 'SouthEast');
    pause
end