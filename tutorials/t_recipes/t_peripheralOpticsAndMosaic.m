function t_peripheralOpticsAndMosaic

    % Examined eccentricities
   
    applyCentralCorrection = ~true;
    
    % Generate optics using the mean Zernike coefficients
    theSubjectIndex = [];  % mean over all subjects
    desiredPupilDiamMM = 3.0;

    eccXrange = [-25 -10 -5 -2 0 2 5 10 25];
    eccYrange =  [-25 -10 -5 0 5 10 25];
    theSubjectIndex = 2;
    %computeOIAndMosaicAcrossEccentricitiesPolans(theSubjectIndex, desiredPupilDiamMM, eccXrange, eccYrange, applyCentralCorrection);
    
    
    whichEye = 'right';
    whichEcc = 'temporal';
    % Eccentricities are specified in visual field space
    temporalEcc = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 10 20];  % Temporal side, which includes the OD at around +(12-15) degs
    nasalEcc = -fliplr(temporalEcc);                                       % Nasal size, which includes the fovea at around -5 deg
    if (strcmp(whichEcc, 'temporal'))
        eccXrange = temporalEcc;
    else
        eccXrange = nasalEcc;
    end
    computeOIAndMosaicAcrossEccentricitiesArtal(desiredPupilDiamMM, eccXrange, whichEye, applyCentralCorrection, sprintf('%s_%s', whichEye, whichEcc));
    
end

function computeOIAndMosaicAcrossEccentricitiesArtal(desiredPupilDiamMM, eccXrange, whichEye, applyCentralCorrection, label)
    
    
    
    % Get a struct with the Polans data
    d = Artal2012Data(applyCentralCorrection);
    
    wavefrontSpatialSamples = 201;
    wavelengthsListToCompute = [450:100:750];
    
    visualizePSFOverThisSpatialSupportArcMin = 5; %13;
    visualizePSTAtThisWavelength = 550;
    

    if (strcmp(whichEye,  'right'))
        examinedSubjects = [...
            1   5  7 16 20; ...
            21 39 42 48 50; ...
            54 57 59 60 62; ...
            64 65 69 70 72; ...
            81 85 90 97 100 ...
        ];
    else
        examinedSubjects = [ ...
            1   5   7  8 16; ...
            20 21  26 29 30; ...
            37 38 39 42 48; ...
            54 62 65 69 70; ...
            79 81 88 97 100];
    end
    
    % plotting coords
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', size(examinedSubjects,2), ...
       'colsNum', numel(eccXrange), ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.01, ...
       'topMargin',      0.01);
   
   

   for subjectGroupIndex = 1:size(examinedSubjects,1)
       
       hFig = figure(subjectGroupIndex+10); clf;
       set(hFig, 'Position', [10 10 4040 1420], 'Color', [1 1 1]);
   
  
   for sIndex = 1:size(examinedSubjects,2)
   for eccXindex = 1:numel(eccXrange)
      
       
       % The eccentricity in degrees
        eccXY = [eccXrange(eccXindex) 0];
        
        % Get cone spacing at this eccentricity
        eccRadiusDegs = sqrt(sum(eccXY.^2,2));
        eccAngleDegs = atan2d(eccXY(2), eccXY(1));
        [coneSpacingInMeters, coneApertureInMeters] = coneSizeReadData('eccentricity', eccRadiusDegs, ...
                                        'angle', eccAngleDegs, ...
                                        'eccentricityUnits', 'deg', ...
                                        'angleUnits', 'deg', ...
                                        'whichEye', whichEye, ...
                                        'useParfor', false);
                                    
        % Generate a 0.2 x 0.2 deg regular hex cone mosaic with ecc-adjusted cone separation and aperture
%         theConeMosaic = coneMosaicHex(13, ...
%             'fovDegs', 0.3, ...
%             'customLambda', coneSpacingInMeters*1e6, ...
%             'customInnerSegmentDiameter', coneApertureInMeters*1e6);


        % Get zCoeffs for this eccentricity
        theSubjectIndex = examinedSubjects(subjectGroupIndex,sIndex);
        [zCoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d,theSubjectIndex, eccXY);
        
        
        % Generate oi at this eccentricity
        theOI = makeCustomOI(zCoeffs, d.measurementPupilDiameMM, d.measurementWavelength, ...
            desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, whichEye);

        % Plot PSF and cone mosaic at this eccentricity
        ax1 = subplot('Position', subplotPosVectors(mod(sIndex-1,5)+1,eccXindex).v);
%         visualizePSF(theOI, visualizePSTAtThisWavelength, visualizePSFOverThisSpatialSupportArcMin, ...
%             'withSuperimposedMosaic', theConeMosaic, ...
%             'contourLevels', [0.1 0.25 0.5 0.75 0.9], ...
%             'includePupilAndInFocusWavelengthInTitle', (eccXindex==1)&&(sIndex == numel(examinedSubjects)), ...
%             'axesHandle', ax1, 'fontSize', 14);
        
        visualizePSF(theOI, visualizePSTAtThisWavelength, visualizePSFOverThisSpatialSupportArcMin, ...
            'contourLevels', [0.1 0.25 0.5 0.75 0.9], ...
            'includePupilAndInFocusWavelengthInTitle', (eccXindex==1)&&(sIndex == numel(examinedSubjects)), ...
            'axesHandle', ax1, 'fontSize', 14);
        
        
        if (any(isnan(zCoeffs)))
            title('NaN values');
        end
        
        if (sIndex < numel(examinedSubjects))
            xlabel('');
            set(gca, 'XTickLabel', {});
        end
        
        if (eccXrange(eccXindex) > min(eccXrange))
            set(gca, 'YTickLabel', {});
        else
            ylabel(sprintf('subject %d', theSubjectIndex));
        end
        drawnow;
        
   end
   end
   NicePlot.exportFigToPDF(sprintf('%s_group%d',label,subjectGroupIndex), hFig, 300);
   close(hFig);
   end
   
end


function computeOIAndMosaicAcrossEccentricitiesPolans(theSubjectIndex, desiredPupilDiamMM, eccXrange, eccYrange, applyCentralCorrection)

    % Get a struct with the Polans data
    d = Polans2015Data(applyCentralCorrection);

    
    wavefrontSpatialSamples = 201;
    wavelengthsListToCompute = [450:100:750];
    
    visualizePSFOverThisSpatialSupportArcMin = 7;
    visualizePSTAtThisWavelength = 550;
    
    % plotting coords
    posVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', numel(eccYrange), ...
       'colsNum', numel(eccXrange), ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.01);
   
   
    if (isempty(theSubjectIndex))
        figNo = 1000;
    else
        figNo = theSubjectIndex;
    end
    
    % Reset figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1440 1180], 'Color', [1 1 1]);
    
    for eccYindex = 1:numel(eccYrange)
    for eccXindex = 1:numel(eccXrange)
    
        % The eccentricity in degrees
        eccXY = [eccXrange(eccXindex) eccYrange(eccYindex)];
        
        % Get cone spacing at this eccentricity
        eccRadiusDegs = sqrt(sum(eccXY.^2,2));
        eccAngleDegs = atan2d(eccXY(2), eccXY(1));
        [coneSpacingInMeters, coneApertureInMeters] = coneSizeReadData('eccentricity', eccRadiusDegs, ...
                                        'angle', eccAngleDegs, ...
                                        'eccentricityUnits', 'deg', ...
                                        'angleUnits', 'deg', ...
                                        'whichEye', d.eye, ...
                                        'useParfor', false);
                                    
        % Generate a 0.2 x 0.2 deg regular hex cone mosaic with ecc-adjusted cone separation and aperture
        theConeMosaic = coneMosaicHex(13, ...
            'fovDegs', 0.3, ...
            'customLambda', coneSpacingInMeters*1e6, ...
            'customInnerSegmentDiameter', coneApertureInMeters*1e6);


        % Get zCoeffs for this eccentricity
        [zCoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d,theSubjectIndex, eccXY);

        % Generate oi at this eccentricity
        theOI = makeCustomOI(zCoeffs, d.measurementPupilDiameMM, d.measurementWavelength, ...
            desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, d.eye);

        % Plot PSF and cone mosaic at this eccentricity
        ax1 = subplot('Position', posVectors(numel(eccYrange)-eccYindex+1,eccXindex).v);
        visualizePSF(theOI, visualizePSTAtThisWavelength, visualizePSFOverThisSpatialSupportArcMin, ...
            'withSuperimposedMosaic', theConeMosaic, ...
            'contourLevels', [0.1 0.25 0.5 0.75 0.9], ...
            'includePupilAndInFocusWavelengthInTitle', (eccXindex==1)&&(eccYindex==1), ...
            'axesHandle', ax1, 'fontSize', 12);
        
        if (eccYrange(eccYindex) > min(eccYrange))
            xlabel('');
            set(gca, 'XTickLabel', {});
        end
        
        if (eccXrange(eccXindex) > min(eccXrange))
            set(gca, 'YTickLabel', {});
        end
        drawnow;
    end
    end
    
    NicePlot.exportFigToPDF(sprintf('PolansSubject%d_pupilDiam%2.1fmm.pdf', figNo, desiredPupilDiamMM), hFig, 300);
    

%     for subjectIndex = 1:d.subjectsNum
%         generateFigure3(d, subjectIndex,nan);
%     end
    
    % Generate figure 3 of the Polans (2015) paper
    
    %generateFigure3(d, theSubjectIndex);
    %generateFigureS3(d, theSubjectIndex)
    
end

function theOI = makeCustomOI(zCoeffs, measPupilDiameterMM, measWavelength, ...
    desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, whichEye)

    showTranslation = false;
    
    [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measPupilDiameterMM, desiredPupilDiamMM, ...
             measWavelength, showTranslation);
        
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(theOTF(:,:,waveIndex));
        theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
    end
       
    umPerDegree = 300;
    theOI = oiCreate('wvf human', desiredPupilDiamMM,[],[], umPerDegree);
    optics = oiGet(theOI,'optics');
    optics = opticsSet(optics, 'otfwave', wavelengthsListToCompute);
    
    % Update optics with new OTF data
    xSfCyclesPerMM = 1000*xSfCyclesDeg / umPerDegree;
    ySfCyclesPerMM = 1000*ySfCyclesDeg / umPerDegree;
    customOptics = opticsSet(optics,'otf data',theOTF);
    customOptics = opticsSet(customOptics, 'otffx',xSfCyclesPerMM);
    customOptics = opticsSet(customOptics,'otffy',ySfCyclesPerMM);
    
    % Update theOI with custom optics
    theOI = oiSet(theOI,'optics', customOptics);
    if (strcmp(whichEye, 'right'))
        if (nearestEccXY(1)<0)
            xEccLabel = sprintf('x:%2.0f^o(N)', nearestEccXY(1));
        elseif(nearestEccXY(1)>0)
            xEccLabel = sprintf('x:%2.0f^o(T)', nearestEccXY(1));
        elseif(nearestEccXY(1)==0)
            xEccLabel = sprintf('x:%2.0f^o', nearestEccXY(1));
        end
    elseif (strcmp(whichEye, 'left'))
        if (nearestEccXY(1)<0)
            xEccLabel = sprintf('x:%2.0f^o(T)', nearestEccXY(1));
        elseif(nearestEccXY(1)>0)
            xEccLabel = sprintf('x:%2.0f^o(N)', nearestEccXY(1));
        elseif(nearestEccXY(1)==0)
            xEccLabel = sprintf('x:%2.0f^o', nearestEccXY(1));
        end
    end
    yEccLabel = sprintf('y:%2.0f^o', nearestEccXY(2));
    theOI = oiSet(theOI,'name', sprintf('%s  %s', xEccLabel, yEccLabel));
    
end

function [theZcoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d,subjectIndex, eccXY)
    dX = d.eccXgrid - eccXY(1);
    dY = d.eccYgrid - eccXY(2);
    [~,indexOfNearestEcc] = min(dX.^2+dY.^2);
    nearestEccXY = [d.eccXgrid(indexOfNearestEcc) d.eccYgrid(indexOfNearestEcc)];
    
    if (strcmp(d.source, 'Polans_2015'))
        if (isempty(subjectIndex))
            % mean over all subjects
            theMeasuredZcoeffs = mean(squeeze(d.zCoeffs(:,indexOfNearestEcc,:)),1);
        else
            theMeasuredZcoeffs = squeeze(d.zCoeffs(subjectIndex,indexOfNearestEcc,:));
        end
    else
        eyeIndex = 1;
        if (isempty(subjectIndex))
            theMeasuredZcoeffs = squeeze(d.zCoeffsMean(eyeIndex,indexOfNearestEcc,:));
        else
            theMeasuredZcoeffs = squeeze(d.zCoeffs(eyeIndex,subjectIndex,indexOfNearestEcc,:));
        end
        
    end
    
    % Place zCoeffs in right bin
    theZcoeffs = zeros(1,21);
    theZcoeffs(d.zCoeffOSAIndices+1) = theMeasuredZcoeffs;
end



function generateFigureS3(d, theSubjectIndex)
    visualizedZCoeffs = {'Z3', 'Z4', 'Z5', 'Z7', 'Z12'};
    allCoeffIndices = d.zCoeffNames(:,1);
    figure(); clf;
    for k = 1:numel(visualizedZCoeffs)
        subplot(2,3,k);
        % Extract data set index for visualized coeffectient
        theDataSetIndex = find(strcmp(allCoeffIndices, visualizedZCoeffs{k}));
        if (~isempty(theDataSetIndex))
            % plot the distribution across the vertical meridian for this coefficient, 
            % across the 10 subjects
            eccIndicesAcrossVerticalMeridian = find(d.eccXgrid == 0);
            eccY = d.eccYgrid(eccIndicesAcrossVerticalMeridian);
            theZcoeffs = squeeze(d.zCoeffs(1:d.subjectsNum,eccIndicesAcrossVerticalMeridian,theDataSetIndex));
            
            if (~isnan(theSubjectIndex))
                % get data across the vertical meridian
                indexOfZeroEcc = find(eccY == 0);
                for eccIndex = 1:numel(eccY)
                    subjectRange = [min(squeeze(theZcoeffs(:,eccIndex))) max(squeeze(theZcoeffs(:,eccIndex)))];
                    subjectMean(eccIndex) = mean(squeeze(theZcoeffs(:,eccIndex)));
                    subjectStd = std(squeeze(theZcoeffs(:,eccIndex)));
                    plot(eccY(eccIndex)*[1 1], subjectMean(eccIndex)+0.5*subjectStd*[-1 1], 'b-');  hold on;
                end
                plot(eccY, subjectMean, 'b-');
                % Bestsubject
                plot(eccY, theZcoeffs(theSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                title(sprintf('%s (%s)', d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));

            else
                % Individual subject
                theSubjectIndex = theSubjectIndex; 
                plot(eccX, theZcoeffs(theSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                hold on
                indexOfZeroEcc = find(eccY == 0);
                if strcmp(d.zCoeffNames{theDataSetIndex,2}, 'defocus')
                    plot(eccY, theZcoeffs(theSubjectIndex,:)-theZcoeffs(theSubjectIndex,indexOfZeroEcc), 'r--', 'LineWidth', 1.5);
                end
                legend({'measured', 'central correction subtracted'});
                title(sprintf('SUBJECT %d %s (%s)', theSubjectIndex, d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));

            end
            
            switch (visualizedZCoeffs{k})
                case 'Z3'
                    set(gca, 'YLim', [-0.5 0.5]);
                case 'Z4'
                    set(gca, 'YLim', [-3 2]);
                case 'Z5'
                    set(gca, 'YLim', [-5 1]);
                case 'Z7'
                    set(gca, 'YLim', [-1 1]);
                case 'Z12'
                    set(gca, 'YLim', [-0.05 0.2]);
            end
            xlabel('<- superior   eccentricity, y (degs)   inferior->');
            ylabel('Zernike aberration (microns)');
        end
    end
end


function generateFigure3(d, theSubjectIndex)
    visualizedZCoeffs = {'Z3', 'Z4', 'Z5', 'Z8', 'Z9', 'Z12'};
    allCoeffIndices = d.zCoeffNames(:,1);
    
    figure(); clf;
    for k = 1:numel(visualizedZCoeffs)
        subplot(2,3,k);
        % Extract data set index for visualized coeffectient
        theDataSetIndex = find(strcmp(allCoeffIndices, visualizedZCoeffs{k}));
        if (~isempty(theDataSetIndex))
            % plot the distribution across the horizontal meridian for this coefficient, 
            % across the 10 subjects
            eccIndicesAcrossHorizontalMeridian = find(d.eccYgrid == 0);
            eccX = d.eccXgrid(eccIndicesAcrossHorizontalMeridian);
            theZcoeffs = squeeze(d.zCoeffs(1:d.subjectsNum,eccIndicesAcrossHorizontalMeridian,theDataSetIndex));
            
            if (~isnan(theSubjectIndex))
                % get data across the horizontal meridian
                for eccIndex = 1:numel(eccX)
                    subjectRange = [min(squeeze(theZcoeffs(:,eccIndex))) max(squeeze(theZcoeffs(:,eccIndex)))];
                    subjectMean(eccIndex) = mean(squeeze(theZcoeffs(:,eccIndex)));
                    subjectStd = std(squeeze(theZcoeffs(:,eccIndex)));
                    plot(eccX(eccIndex)*[1 1], subjectMean(eccIndex)+0.5*subjectStd*[-1 1], 'b-');  hold on;
                end
                plot(eccX, subjectMean, 'b-');
                % Best subject
                plot(eccX, theZcoeffs(theSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                title(sprintf('%s (%s)', d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));

            else
                % Individual subject
                theSubjectIndex = theSubjectIndex; 
                plot(eccX, theZcoeffs(theSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                hold on
                indexOfZeroEcc = find(eccX == 0);
                if strcmp(d.zCoeffNames{theDataSetIndex,2}, 'defocus')
                    plot(eccX, theZcoeffs(theSubjectIndex,:)-theZcoeffs(theSubjectIndex,indexOfZeroEcc), 'r--', 'LineWidth', 1.5);
                end
                legend({'measured', 'central correction subtracted'});
                title(sprintf('SUBJECT %d %s (%s)', theSubjectIndex, d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));
            end
            % add box around the optic disk (10 to 18 degs)
            b.x = -14 +[-4 -4 4 4 4];
            b.y = [-1.5 1.5 1.5 -1.5 -1.5];
            patch(b.x, b.y, [0.5 0.5 0.5], 'FaceAlpha', 0.5);
            
            set(gca, 'YLim', [-1.5 1.5]);
            xlabel('<- nasal   eccentricity, x (degs)   temporal->');
            ylabel('Zernike aberration (microns)');
        end
    end
end

function d = Artal2012Data(applyCentralCorrection)
    % Load the matfile
    allData = rawDataReadData('zCoefsJaekenArtal2012', ...
                    'datatype', 'isetbiomatfileonpath');
                
    allData = allData.data;
    
    d.source = 'Artal_2012';
    
    d.subjectsNum = size(allData,1);
    
    % Retrieve the Y-eccentricity grid          
    subjectIndex = 1; entryIndex = 4;
    d.eccXgrid = allData(subjectIndex,entryIndex:end);
    d.eccYgrid = 0*d.eccXgrid;
    
    % Truncate headers and reshape data
    allData = allData(2:end, 4:end);
    
    % 2*(130*2)*15 x 81 eccentricities
    d.subjectsNum  = 130;
    d.eyes = {'right', 'left'};
    
    
    % Add zcoeff OSA indices and names
    d.zCoeffOSAIndices = 0:14;
    d.zCoeffNames = {...
          'Z0',  '?'; ...
          'Z1',  '?'; ...
          'Z2',  '?'; ...
          'Z3',  'oblique astigmatism'; ...
          'Z4',  'defocus'; ...
          'Z5',  'vertical astigmatism'; ...
          'Z6' , 'vertical trefoil'; ...
          'Z7',  'vertical coma'; ...
          'Z8',  'horizontal coma'; ...
          'Z9',  'oblique trefoil'; ...
          'Z10', 'oblique quadrafoil'; ...
          'Z11', 'oblique secondary astigmatism'; ...
          'Z12', 'primary spherical'; ...
          'Z13', 'vertical secondary astigmatism'; ...
          'Z14', 'vertical quadrafoil' ...
        };
    
    % Arrange zCoeffs as [eye subject eccentricity ZcoffOrder]
    allData = reshape(allData, [numel(d.zCoeffOSAIndices) d.subjectsNum numel(d.eyes) numel(d.eccXgrid)]);
    zCoeffDimIndex = 1;
    subjectDimIndex = 2;
    eyeDimIndex = 3;
    eccDimIndex = 4;
    d.zCoeffs = permute(allData, [eyeDimIndex subjectDimIndex eccDimIndex zCoeffDimIndex]);
    d.zCoeffsMean = squeeze(mean(d.zCoeffs, 2, 'omitnan'));
    
    % Find subjects whose coeffs are nan and remove them
    tmpZcoeffs = [];
    keptEccIndices = find(abs(d.eccXgrid)<=1);
    
    for subjectIndex = 1:d.subjectsNum
        left = squeeze(d.zCoeffs(2,subjectIndex,keptEccIndices,:));
        if (any(isnan(left(:))))
            fprintf('Removing left eye of suject %d\n', subjectIndex);
        end
        right = squeeze(d.zCoeffs(1,subjectIndex,keptEccIndices,:));
        if (any(isnan(right(:))))
            fprintf('Removing right eye of suject %d\n', subjectIndex);
        end
    end

    % Using a 4 mm pupil
    d.measurementPupilDiameMM = 4;
    
    % And 780 nm light
    d.measurementWavelength = 550;
    
    if (applyCentralCorrection)
        d.zCoeffs = removeCentralRefractionError(d);
    end
    
    function zCoeffs = removeCentralRefractionError(d)
        zCoeffs = d.zCoeffs;
        % Find the zero ecc index
        zeroEccIndex = find(d.eccXgrid == 0);
        % Find the z-coeff index corresponding to the defocus term
        theDefocusCoeffIndex = find(strcmp(d.zCoeffNames(:,2), 'defocus'));
        % Central corrections for all subjects
        allSubjectsCentralRefractionCorrection = zCoeffs(:, :, zeroEccIndex, theDefocusCoeffIndex);
        % Subtract central correction for all subjects and all eccentricities
        zCoeffs(:,:,:,theDefocusCoeffIndex) = bsxfun(@minus, ...
            zCoeffs(:,:,:,theDefocusCoeffIndex), allSubjectsCentralRefractionCorrection);
    end

end

function d = Polans2015Data(applyCentralCorrection)
% Load the raw data from Polans et al (2015), "Wide-field optical model
%     of the human eye with asymmetrically tilted and decentered lens 
%     that reproduces measured ocular aberrations", Optica, 2(2), 2015,
%     pp.124-134
%
    
    % Load the matfile
    allData = rawDataReadData('zCoefsPolans2015', ...
                        'datatype', 'isetbiomatfileonpath');
    allData = allData.data;
    
    d.source = 'Polans_2015';
    d.subjectsNum = size(allData,1);
    % Retrieve the Y-eccentricity grid          
    subjectIndex = 1; entryIndex = 1;
    d.eccYgrid = allData(subjectIndex,:,entryIndex);
    
    
    % Retrieve the X-eccentricity grid
    entryIndex = 2;
    d.eccXgrid = allData(subjectIndex,:,entryIndex);
    
    % Retrieve the z-coeffs
    entryIndices = 3:20;                
    d.zCoeffs = allData(:,:,entryIndices);
    if (any(isnan(d.zCoeffs(:))))
        fprintf(2,'Found z-coeff with NaN values');
    end

    
    % Add zcoeff OSA indices and names
    d.zCoeffOSAIndices = entryIndices;
    d.zCoeffNames = {...
          'Z3',  'oblique astigmatism'; ...
          'Z4',  'defocus'; ...
          'Z5',  'vertical astigmatism'; ...
          'Z6' , 'vertical trefoil'; ...
          'Z7',  'vertical coma'; ...
          'Z8',  'horizontal coma'; ...
          'Z9',  'oblique trefoil'; ...
          'Z10', 'oblique quadrafoil'; ...
          'Z11', 'oblique secondary astigmatism'; ...
          'Z12', 'primary spherical'; ...
          'Z13', 'vertical secondary astigmatism'; ...
          'Z14', 'vertical quadrafoil'; ...
          'Z15', '??'; ...
          'Z16', '??'; ...
          'Z17', '??'; ...
          'Z18', '??'; ...
          'Z19', '??'; ...
          'Z20', '??' ...
        };
    
    % Measurement params
   
    % The Polans 2015 measurements were all done in the right eye
    d.eye = 'right';
    
    % Using a 4 mm pupil
    d.measurementPupilDiameMM = 4;
    
    % And 550 nm light
    d.measurementWavelength = 550;
    
    if (applyCentralCorrection)
        d.zCoeffs = removeCentralRefractionError(d);
    end
    
    function zCoeffs = removeCentralRefractionError(d)
        zCoeffs = d.zCoeffs;
        % Find the zero ecc index
        zeroEccIndex = find(d.eccXgrid == 0 & d.eccYgrid == 0);
        % Find the z-coeff index corresponding to the defocus term
        theDefocusCoeffIndex = find(strcmp(d.zCoeffNames(:,2), 'defocus'));
        % Central corrections for all subjects
        allSubjectsCentralRefractionCorrection = zCoeffs(:, zeroEccIndex, theDefocusCoeffIndex);
        % Subtract central correction for all subjects and all eccentricities
        zCoeffs(:,:,theDefocusCoeffIndex) = bsxfun(@minus, ...
            zCoeffs(:,:,theDefocusCoeffIndex), allSubjectsCentralRefractionCorrection);
    end

end

