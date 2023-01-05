function photonCatchComputations()

    assembleISETBioData = ~true;
    if (assembleISETBioData)
        % You need to have ISETBio in your path to run this
        assembleISETBioDataAcrossEccentricities();
    end

    % You need to have ISETBioData.mat in your path to run this
    analyzePhotonCatchFactors();
end

function assembleISETBioDataAcrossEccentricities()
    mappedRFsDir = '/Volumes/SSDdisk/MATLAB/toolboxes/isetbio/isettools/ganglioncells/JohannesAnalysesDataOLD';
    ZernikeDataBase = 'Polans2015';
    subjectRankOrder = 6;

    % Examined eccentricities
    eccX = -8:1:8;
    eccY = -6:1:6;
    [eccXGrid, eccYGrid] = meshgrid(eccX, eccY);
    eccDegsGrid = [eccXGrid(:) eccYGrid(:)];

    dataDict = containers.Map();

    % Eccentricity data
    for iEcc = 1:size(eccDegsGrid,1)
        dataLabel = sprintf('eccXY = %2.2f,%2.2f', eccDegsGrid(iEcc,1), eccDegsGrid(iEcc,2));
        dataDict(dataLabel) = retrieveISETBioComponents(mappedRFsDir, ZernikeDataBase, subjectRankOrder, eccDegsGrid(iEcc,:));
    end

    % Meta data
    dataDict('metaData') = struct(...
        'LCONE_ID', cMosaic.LCONE_ID, ...
        'MCONE_ID', cMosaic.MCONE_ID, ...
        'SCONE_ID', cMosaic.SCONE_ID, ...
        'cMap', cat(1,brewermap(512,'*YlGnBu'), brewermap(512,'YlOrRd')) ...
        );

    % Save assembled data
    save('ISETBioData.mat', 'eccDegsGrid',  'dataDict');
end

function dTarget = retrieveISETBioComponents(mappedRFsDir, ZernikeDataBase, subjectRankOrder, targetEccDegs)
    % Retrieve data
    innerSegmentDiameterDegsAtTargetEccentricity = ...
        innerSegmentDiameterDegsAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, targetEccDegs);

    [outerSegmentAxialLengthMicronsAtTargetEccentricity, absorptanceSpectraAtTargetEccentricity] = ...
        outerSegmentLengthMicronsAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, targetEccDegs);

    macularPigmentTransmittanceAtTargetEccentricity = ...
        macularPigmentTransmittanceAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, targetEccDegs);

    meanConesNumInRFcenterAtTargetEccentricity = ...
        meanConesNumInMRGCRFcenterAtTargetEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, targetEccDegs);

    % Target ecc data struct
    fprintf('Assembling data for %f %f degs\n', targetEccDegs(1), targetEccDegs(2));
    dTarget = struct(...
        'innerSegmentDiameterDegs', innerSegmentDiameterDegsAtTargetEccentricity , ...
        'outerSegmentLengthMicrons', outerSegmentAxialLengthMicronsAtTargetEccentricity, ...
        'absorptanceSpectra', absorptanceSpectraAtTargetEccentricity, ...
        'macularPigmentTransmittance', macularPigmentTransmittanceAtTargetEccentricity, ...
        'meanConesNumInRFcenter', meanConesNumInRFcenterAtTargetEccentricity);
end

function meanConesNumInRFcenterAtTargetEccentricity = ...
        meanConesNumInMRGCRFcenterAtTargetEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, mosaicEccDegs)

    % Load the computed components data for the fovea
    fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
    fName = fullfile(mappedRFsDir, fName);
    load(fName, 'theMidgetRGCmosaic');

    % Mean number of cones in RF center at the center of the mosaic
    rgcDistancesFromMosaicCenter = sum((bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, mosaicEccDegs)).^2,2);
    [~,centerMostRGCIndices] = sort(rgcDistancesFromMosaicCenter, 'ascend');
    centerMostRGCsNum = 10;
    centerMostRGCIndices = centerMostRGCIndices(1:centerMostRGCsNum);
    meanConesNumInRFcenterAtTargetEccentricity = mean(sum(full(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix(:, centerMostRGCIndices)),1));
end


function macularPigmentTransmittanceAtTargetEccentricity = ...
    macularPigmentTransmittanceAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, mosaicEccDegs)

    % Load the computed components data for the fovea
    fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
    fName = fullfile(mappedRFsDir, fName);
    load(fName, 'theMidgetRGCmosaic');

    % Retrieve the macular pigment object
    theMacular = theMidgetRGCmosaic.inputConeMosaic.macular;

    % Compute transmittance at target eccentricity
    macularPigmentDensityAtMosaicEcc = theMacular.eccDensity([], 'eccDegs2', sum(mosaicEccDegs.^2,2));
    macularPigmentTransmittanceAtTargetEccentricity = 10.^(-macularPigmentDensityAtMosaicEcc * theMacular.unitDensity');
end


function [coneOuterSegmentAxialLengthMicronsAtTargetEccentricity, absorptanceSpectraAtTargetEccentricity] = ...
    outerSegmentLengthMicronsAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, mosaicEccDegs)

    % Load the computed components data for the fovea
    fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, 0,0);
    fName = fullfile(mappedRFsDir, fName);
    load(fName, 'theMidgetRGCmosaic');

    % Retrieve the photopigment object
    thePhotoPigment = theMidgetRGCmosaic.inputConeMosaic.pigment;

    % Specific density of foveal L-cones is 0.013 +/- 0.002 per micron for the L-cones (Bowmaker et al., 1978)
    LconeType = cMosaic.LCONE_ID;
    LConeSpecificDensity = 0.013;
    
    % Compute the axial length of the foveal outer segment
    fovealLConeAxialOpticalDensity = thePhotoPigment.opticalDensity(LconeType);
    fovealConeOuterSegmentAxialLengthMicrons = fovealLConeAxialOpticalDensity / LConeSpecificDensity;

    % Comute the cone specific densities
    for coneTypeIndex = cMosaic.LCONE_ID: cMosaic.SCONE_ID
        coneSpecificDensity(coneTypeIndex) = thePhotoPigment.opticalDensity(coneTypeIndex)/fovealConeOuterSegmentAxialLengthMicrons;
    end

    % Load the computed components data at the target ecc
    fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
    fName = fullfile(mappedRFsDir, fName);
    load(fName, 'theMidgetRGCmosaic');

     % Retrieve the photopigment object
    thePhotoPigment = theMidgetRGCmosaic.inputConeMosaic.pigment;

    % Compute the target cone outer segment axial length
    coneOuterSegmentAxialLengthMicronsAtTargetEccentricity = fovealConeOuterSegmentAxialLengthMicrons * ...
        outerSegmentLengthVariationWithEccentricity(mosaicEccDegs);

    % Compute the  L-, M-, and S-cone absorptance spectra at the target eccentricity
    for coneTypeIndex = cMosaic.LCONE_ID: cMosaic.SCONE_ID
        axialOpticalDensityAtTargetEccentricity = coneSpecificDensity(coneTypeIndex) * coneOuterSegmentAxialLengthMicronsAtTargetEccentricity;
        absorbanceSpectrum = thePhotoPigment.absorbance(:,coneTypeIndex);
        absorptanceSpectraAtTargetEccentricity(coneTypeIndex,:) = 1 - 10 .^ (-absorbanceSpectrum*axialOpticalDensityAtTargetEccentricity);
    end
end


function innerSegmentDiameterDegsAtTargetEccentricity = innerSegmentDiameterDegsAtEccentricity(mappedRFsDir, ZernikeDataBase, subjectRankOrder, mosaicEccDegs)
    % Load the computed components data
    fName = sprintf('mosaicAnd%s_Subject%d_optics_EccXY_%2.2f_%2.2f.mat', ...
        ZernikeDataBase, subjectRankOrder, mosaicEccDegs(1), mosaicEccDegs(2));
    fName = fullfile(mappedRFsDir, fName);
    load(fName, 'theMidgetRGCmosaic');

    % Sort cones according to their distance from the mosaic center
    d = bsxfun(@minus, theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs, mosaicEccDegs);
    r = sum(d.^2,2);
    [~,sortedConeIndices] = sort(r, 'ascend');

    % Take the mean diameter across 6 cones at the target eccentricity
    innerSegmentDiameterDegsAtTargetEccentricity = mean(theMidgetRGCmosaic.inputConeMosaic.coneApertureDiametersDegs(sortedConeIndices(1:6)));
end


function relativeWithRespectToFoveaOuterSegmentLength = outerSegmentLengthVariationWithEccentricity(eccDegs)
% Scanned data (eccentricity, osLength in microns) from Figure 1 (right panel)
% Banks, Sekuler and Anderson (1991). Peripher spatial vision: limits
% imposed by optics, photoreceptors and receptor pooling
    s = [ ...
        0.00  47.81;
        1.82  26.16;
        4.86  21.2;
        9.86  21.20;
        19.78  21.2;
        39.90  13.22;
    ];
  switch (size(eccDegs,2))
      case 1 
          % do nothing
      case 2
          eccDegs = sqrt(sum(eccDegs.^2,2));
      otherwise
          error('eccDegs must be either [nx1] or [nx2]')
  end
  scannedData.eccDegsRaw = s(:,1);
  scannedData.lengthMicronsRaw = s(:,2);
  interpolationMethod = 'pchip';
  osLengthMicrons = interp1(scannedData.eccDegsRaw, scannedData.lengthMicronsRaw, eccDegs, interpolationMethod);
  osLengthAtZeroEcc = s(1,2);
  relativeWithRespectToFoveaOuterSegmentLength = osLengthMicrons / osLengthAtZeroEcc;
end

