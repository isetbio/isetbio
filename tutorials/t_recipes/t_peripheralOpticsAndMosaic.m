function t_peripheralOpticsAndMosaic

    % Examined eccentricities
    eccXrange = -40:20:40;
    eccYrange =  sort(-10:5:10, 'descend');
 
    % Generate optics using the mean Zernike coefficients
    theSubjectIndex = [];  % mean over all subjects
    computeOIAndMosaicAcrossEccentricities(theSubjectIndex, desiredPupilDiamMM, eccXrange, eccYrange);
    
    % Generate optics using individual subject Zernike coefficients
    for theSubjectIndex = 1:10
        desiredPupilDiamMM = 3.0;
        computeOIAndMosaicAcrossEccentricities(theSubjectIndex, desiredPupilDiamMM, eccXrange, eccYrange);
    end
end

function computeOIAndMosaicAcrossEccentricities(theSubjectIndex, desiredPupilDiamMM, eccXrange, eccYrange)

    % Get a struct with the Polans data
    applyCentralCorrection = true;
    d = PolansData(applyCentralCorrection);
     
    wavefrontSpatialSamples = 201;
    wavelengthsListToCompute = 450:100:750;
    
    visualizePSFOverThisSpatialSupportArcMin = 8;
    visualizePSTAtThisWavelength = d.measurementWavelength;
    
    % plotting coords
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', numel(eccYrange), ...
       'colsNum', numel(eccXrange), ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.03, ...
       'topMargin',      0.03);
   
   
    if (isempty(theSubjectIndex))
        figNo = 1000;
    else
        figNo = theSubjectIndex;
    end
    
    % Reset figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1700 1500], 'Color', [1 1 1]);
    
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
            'fovDegs', 0.2, ...
            'customLambda', coneSpacingInMeters*1e6, ...
            'customInnerSegmentDiameter', coneApertureInMeters*1e6);


        % Get zCoeffs for this eccentricity
        [zCoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d,theSubjectIndex, eccXY);

        % Generate oi at this eccentricity
        theOI = makeCustomOI(zCoeffs, d.measurementPupilDiameMM, d.measurementWavelength, ...
            desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, d.eye);

        % Plot PSF and cone mosaic at this eccentricity
        ax1 = subplot('Position', subplotPosVectors(eccYindex,eccXindex).v);
        visualizePSF(theOI, visualizePSTAtThisWavelength, visualizePSFOverThisSpatialSupportArcMin, ...
            'withSuperimposedMosaic', theConeMosaic, ...
            'contourLevels', [0.1 0.25 0.5 0.75 0.9], ...
            'axesHandle', ax1, 'fontSize', 14);
        
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
    
    NicePlot.exportFigToPDF(sprintf('Subject%d_pupilDiam%2.1fmm.pdf', figNo, desiredPupilDiamMM), hFig, 300);
    

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
            xEccLabel = sprintf('%2.0f^o (N)', nearestEccXY(1));
        elseif(nearestEccXY(1)>0)
            xEccLabel = sprintf('%2.0f^o (T)', nearestEccXY(1));
        elseif(nearestEccXY(1)==0)
            xEccLabel = sprintf('%2.0f^o', nearestEccXY(1));
        end
    elseif (strcmp(whichEye, 'left'))
        if (nearestEccXY(1)<0)
            xEccLabel = sprintf('%2.0f^o (T)', nearestEccXY(1));
        elseif(nearestEccXY(1)>0)
            xEccLabel = sprintf('%2.0f^o (N)', nearestEccXY(1));
        elseif(nearestEccXY(1)==0)
            xEccLabel = sprintf('%2.0f^o', nearestEccXY(1));
        end
        theOI = oiSet(theOI,'name', sprintf('%2.0f^o,%2.0f^o', nearestEccXY(1), nearestEccXY(2)));
    end
    theOI = oiSet(theOI,'name', sprintf('%s,%2.0f^o', xEccLabel, nearestEccXY(2)));
    
end

function [theZcoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d,subjectIndex, eccXY)
    dX = d.eccXgrid - eccXY(1);
    dY = d.eccYgrid - eccXY(2);
    [~,indexOfNearestEcc] = min(dX.^2+dY.^2);
    nearestEccXY = [d.eccXgrid(indexOfNearestEcc) d.eccYgrid(indexOfNearestEcc)];
    if (isempty(subjectIndex))
        % mean over all subjects
        theZcoeffs = mean(squeeze(d.zCoeffs(:,indexOfNearestEcc,:)),1);
    else
        theZcoeffs = squeeze(d.zCoeffs(subjectIndex,indexOfNearestEcc,:));
    end
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


function d = PolansData(applyCentralCorrection)
% Load the raw data from Polans et al (2015), "Wide-field optical model
%     of the human eye with asymmetrically tilted and decentered lens 
%     that reproduces measured ocular aberrations", Optica, 2(2), 2015,
%     pp.124-134
%
    
    % Load the matfile
    allData = rawDataReadData('zCoefsPolans2015', ...
                        'datatype', 'isetbiomatfileonpath');
    allData = allData.data;
    
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
        theDataSetIndex = find(strcmp(d.zCoeffNames(:,2), 'defocus'));
        % Central corrections for all subjects
        allSubjectsCentralRefractionCorrection = zCoeffs(:, zeroEccIndex, theDataSetIndex);
        % Subtract central correction for all subjects and all eccentricities
        zCoeffs(:,:,theDataSetIndex) = bsxfun(@minus, ...
            zCoeffs(:,:,theDataSetIndex), allSubjectsCentralRefractionCorrection);
    end

end

