function t_peripheralOpticsAndMosaic

    %% Initialize
    ieInit;

    % Get a struct with the Polans data
    d = PolansData();
    
    bestSubjectIndex = 7;
    worseSubjectIndex = 9;
    
    computeCentralRefractionErrorForAllSubjects(d);
    
    meridian = 'horizontal';
    applyCentralCorrection = true
    
    
    [zCoeffs, eccDegs] = zCoeffsForSubjectAndMeridian(d,bestSubjectIndex, meridian);
    
    targetWavelength = 550;
    targetEccX = -50; 
    targetPupilDiamMM = 2.0;
    theOI = makeCustomOI(zCoeffs, eccDegs, d.measurementPupilDiameMM, d.measurementWavelength, targetEccX, targetPupilDiamMM)

    
    spatialSupportArcMin = 5;
    spatialFrequencySupportCyclePerDeg = 60;
    visualizeOptics(theOI, targetWavelength,spatialSupportArcMin, spatialFrequencySupportCyclePerDeg)

%     for subjectIndex = 1:d.subjectsNum
%         generateFigure3(d, subjectIndex,nan);
%     end
    
    % Generate figure 3 of the Polans (2015) paper
    
    %generateFigure3(d, bestSubjectIndex, worseSubjectIndex);
    %generateFigureS3(d, bestSubjectIndex, worseSubjectIndex)
    
end

function theOI = makeCustomOI(zCoeffsAllEcc, eccDegs, measPupilDiameterMM, measWavelength, targetEccX, targetPupilDiamMM)

    [~,eccIndex] = min(abs(eccDegs-targetEccX));
    zCoeffs = squeeze(zCoeffsAllEcc(eccIndex,:))
    eccDegs(eccIndex)
    
    wavelengthsListToCompute = 400:50:750;
    wavefrontSpatialSamples = 201;
    showTranslation = false;
    
    [thePSF, theOTF, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measPupilDiameterMM, targetPupilDiamMM, ...
             measWavelength, showTranslation);
        
    for waveIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = squeeze(theOTF(:,:,waveIndex));
        theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
    end
       
    umPerDegree = 300;
    theOI = oiCreate('wvf human', targetPupilDiamMM,[],[], umPerDegree);
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
    theOI = oiSet(theOI,'name', sprintf('human-WVF at %2.1f degs ecc', targetEccX));
    
end


function [theZcoeffs, eccDegs] = zCoeffsForSubjectAndMeridian(d, subjectIndex, meridian)
    if (strcmp(meridian, 'horizontal'))
        eccIndicesAcrossMeridian = find(d.eccYgrid == 0);
        eccDegs = d.eccXgrid(eccIndicesAcrossMeridian);
    else
        eccIndicesAcrossMeridian = find(d.eccXgrid == 0);
        eccDegs = d.eccYgrid(eccIndicesAcrossMeridian);
    end
    theZcoeffs = squeeze(d.zCoeffs(subjectIndex,eccIndicesAcrossMeridian,:));
end

function computeCentralRefractionErrorForAllSubjects(d)
    zeroEccIndex = find(d.eccXgrid == 0 & d.eccYgrid == 0);
    theDataSetIndex = find(strcmp(d.zCoeffNames(:,2), 'defocus'));
    centralRefractionCorrections = d.zCoeffs(:, zeroEccIndex, theDataSetIndex);
    for k = 1:d.subjectsNum
        fprintf('Subject %d has central refraction correction (defocus) of : %f\n', ...
            k, centralRefractionCorrections(k));
    end
end

function generateFigureS3(d, bestSubjectIndex, worseSubjectIndex)
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
            
            if (~isnan(bestSubjectIndex) && ~isnan(worseSubjectIndex))
                % get data across the vertical meridian
                indexOfZeroEcc = find(eccY == 0);
                for eccIndex = 1:numel(eccY)
                    subjectRange = [min(squeeze(theZcoeffs(:,eccIndex))) max(squeeze(theZcoeffs(:,eccIndex)))];
                    subjectMean(eccIndex) = mean(squeeze(theZcoeffs(:,eccIndex)));
                    subjectStd = std(squeeze(theZcoeffs(:,eccIndex)));
                    plot(eccY(eccIndex)*[1 1], subjectMean(eccIndex)+0.5*subjectStd*[-1 1], 'b-');  hold on;
                end
                plot(eccY, subjectMean, 'b-');
                % Best and worse subject
                plot(eccY, theZcoeffs(bestSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                plot(eccY, theZcoeffs(worseSubjectIndex,:), 'k-', 'LineWidth', 1.5);
                title(sprintf('%s (%s)', d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));

            else
                % Individual subject
                theSubjectIndex = bestSubjectIndex; 
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


function generateFigure3(d, bestSubjectIndex, worseSubjectIndex)
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
            
            if (~isnan(bestSubjectIndex) && ~isnan(worseSubjectIndex))
                % get data across the horizontal meridian
                for eccIndex = 1:numel(eccX)
                    subjectRange = [min(squeeze(theZcoeffs(:,eccIndex))) max(squeeze(theZcoeffs(:,eccIndex)))];
                    subjectMean(eccIndex) = mean(squeeze(theZcoeffs(:,eccIndex)));
                    subjectStd = std(squeeze(theZcoeffs(:,eccIndex)));
                    plot(eccX(eccIndex)*[1 1], subjectMean(eccIndex)+0.5*subjectStd*[-1 1], 'b-');  hold on;
                end
                plot(eccX, subjectMean, 'b-');
                % Best and worse subject
                plot(eccX, theZcoeffs(bestSubjectIndex,:), 'r-', 'LineWidth', 1.5);
                plot(eccX, theZcoeffs(worseSubjectIndex,:), 'k-', 'LineWidth', 1.5);
                title(sprintf('%s (%s)', d.zCoeffNames{theDataSetIndex,2}, d.zCoeffNames{theDataSetIndex,1}));

            else
                % Individual subject
                theSubjectIndex = bestSubjectIndex; 
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


function d = PolansData()
% Load the raw data from Polans et al (2015), "Wide-field optical model
%     of the human eye with asymmetrically tilted and decentered lens 
%     that reproduces measured ocular aberrations", Optica, 2(2), 2015,
%     pp.124-134
%
    
    % Load the matfile
    allData = rawDataReadData('zCoefsPolans2015', ...
                        'datatype', 'isetbiomatfileonpath');
    allData = allData.data;

    horizontalEccSamplesNum = 81;
    verticalEccSamplesNum = 11;
    
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
    
end

