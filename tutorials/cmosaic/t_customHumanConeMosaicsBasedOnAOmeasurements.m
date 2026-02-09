function t_customHumanConeMosaicsBasedOnAOmeasurements()
% Gnerate ISETBio models of cone mosaics & optics based on AO data from the Sabesan lab
%
% Description:
%   Script used to generate ISETBio models of cone mosaics with corresponding
%   optics based on AO data from the Sabesan lab
%

% History:
%    07/25/25  NPC  ISETBIO Team, Copyright 2025   Wrote it.

    % Specify a source file where the Sabesan data live
    theDataFileName = fullfile(isetbioRootPath, 'isettools/data/cones','RamSabesanData/AO001R_Nasal_v1.csv');

    % Import data
    allMosaicsData = importAllSubjectData(theDataFileName);

    % Generate cone mosaic and optics models for all data
    pupilDiameterMM = []; % generate optics for the measured pupil
    %pupilDiameterMM = 3.0; % generate optics for 3.0 mm pupil
    [theConeMosaics, theHumanOptics, visualizedFOVdegs] = generateModels(allMosaicsData, pupilDiameterMM);

    % Visualize all generated cone mosaics and optics models
    psfVisualizationWavelength = 550;
    visualizeData(theConeMosaics, theHumanOptics, visualizedFOVdegs, psfVisualizationWavelength);
end

function [theConeMosaics,  theHumanOptics, visualizedFOVdegs] = generateModels(allMosaicsData, pupilDiameterMM)

    mosaicsNum = round(size(allMosaicsData,2)/8);
    theConeMosaics = cell(1, mosaicsNum);
    theHumanOptics = cell(1, mosaicsNum);
    fovDegs = zeros(1, mosaicsNum);

    ZernikesCoeffsToInclude = [];
    wavelengthsListToCompute = 400:10:900;

    for iMosaic = 1:mosaicsNum
        [coneMosaicDataTmp, coneMosaicMetaDataTmp, ZernikeDataTmp] = ...
            parseSingleMosaicData(allMosaicsData, iMosaic);

        % Generate an ISETBio model of the cone mosaic based on the imported cone mosaic data
        theConeMosaics{iMosaic} = generateISETBioConeMosaic(...
            coneMosaicDataTmp, coneMosaicMetaDataTmp);

        % Generate an ISETBio model of the optics based on the imported Zernike coefficients
        theHumanOptics{iMosaic} = generateISETBioOptics(...
            ZernikeDataTmp, coneMosaicMetaDataTmp, ...
            pupilDiameterMM, wavelengthsListToCompute, ZernikesCoeffsToInclude);

        fovDegs(iMosaic) = max(coneMosaicMetaDataTmp.fovMM(:)) * 1e3 / coneMosaicMetaDataTmp.retinalMagnificationMicronsPerDegree;
    end

    visualizedFOVdegs = 1.05*max(fovDegs);
end

function visualizeData(theConeMosaics, theHumanOptics, visualizedFOVdegs, psfVisualizationWavelength)
    
    psfVisualizationFovDegs = 0.15;

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1530 1420]);

    sv = NicePlot.getSubPlotPosVectors(...
        'colsNum', numel(theConeMosaics), ...
        'rowsNum', 2, ...
        'widthMargin', 0.03, ...
        'leftMargin', 0.03, ...
        'rightMargin', 0.01, ...
        'heightMargin', 0.07, ...
        'topMargin', 0.02, ...
        'bottomMargin', 0.03);
        
    for iMosaic = 1:numel(theConeMosaics)

        % Retrieve the cone mosaic model
        theSubjectConeMosaic = theConeMosaics{iMosaic};

        % Retrieve the optics model
        theSubjectOI = theHumanOptics{iMosaic};

        % Get the optics data
        optics = oiGet(theSubjectOI, 'optics');

        % Get the PSF data at the visualization wavelength
        psf = opticsGet(optics,'psf data', psfVisualizationWavelength, 'um'); 

        % The PSF amplitude
        thePSFData.data = psf.psf;

        % The PSF spatial support
        focalLength = opticsGet(optics, 'focal length');
        micronsPerDegree = focalLength*tand(1)*1e6;
        thePSFData.supportXdegs = squeeze(psf.xy(1,:,1))/micronsPerDegree;
        thePSFData.supportYdegs = thePSFData.supportXdegs;

        ax = subplot('Position', sv(1,iMosaic).v);
        xo = theSubjectConeMosaic.eccentricityDegs(1);
        yo = theSubjectConeMosaic.eccentricityDegs(2);
        theSubjectConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'noYLabel', (iMosaic>1), ...
            'conesEdgeAlpha', 0.0, ...
            'fontSize', 20, ...
            'domainVisualizationTicks', struct('x', -20:0.5:20, 'y', -20:0.5:20), ...
            'domainVisualizationLimits', [xo-0.5*visualizedFOVdegs xo+0.5*visualizedFOVdegs yo-0.5*visualizedFOVdegs yo+0.5*visualizedFOVdegs]);
        
        
        ax = subplot('Position', sv(2,iMosaic).v);
        xo = theSubjectConeMosaic.eccentricityDegs(1);
        yo = theSubjectConeMosaic.eccentricityDegs(2);
        
        theSubjectConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax, ...
            'conesEdgeAlpha', 0.5, ...
            'conesAlpha', 0.5, ...
            'visualizedConeApertureThetaSamples', 32, ...
            'noYLabel', (iMosaic>1), ...
            'noXLabel', true, ...
            'fontSize', 20, ...
            'withSuperimposedPSF', thePSFData, ...
            'plotTitle', sprintf('measured pupil, %2.0f nm', psfVisualizationWavelength), ...
            'domainVisualizationTicks', struct('x', -20:0.05:20, 'y', -20:0.05:20), ...
            'domainVisualizationLimits', [xo-0.5*psfVisualizationFovDegs xo+0.5*psfVisualizationFovDegs  yo-0.5*psfVisualizationFovDegs yo+0.5*psfVisualizationFovDegs]);
 
        drawnow;
    end % iMosaic

    NicePlot.exportFigToPDF('newformat.pdf', hFig, 300);
end


function theOptics = generateISETBioOptics(theZernikeData, theConeMosaicMetaData, ...
    pupilDiamMM, wavelengthsListToCompute, ZernikesCoeffsToInclude)

    assert(~isempty(strfind(lower(theZernikeData.coeffNames{1}), 'piston')), 'first coeff must be the piston');

    if (~isempty(ZernikesCoeffsToInclude))
        zCoeffs = theZernikeData.coeffValues(1:ZernikesCoeffsToInclude);
    else
        zCoeffs = theZernikeData.coeffValues;
    end

    micronsPerDegree = theConeMosaicMetaData.retinalMagnificationMicronsPerDegree;
    wavefrontSpatialSamples = 301;
   
    measurementPupilDiameterMM = theConeMosaicMetaData.ZernikeCoeffPupilDiameterMM;
    if (isempty(pupilDiamMM))
        pupilDiamMM = measurementPupilDiameterMM;
    end

    inFocusWavelength = 550;
    flipPSFUpsideDown = true;
    upsampleFactor = [];
    noLCA = false;

    % Compute wavefront map
    [~, ~, ~,~, ~, ~, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measurementPupilDiameterMM, ...
             pupilDiamMM, inFocusWavelength, false, ...
             'doNotZeroCenterPSF', false, ...
             'micronsPerDegree', micronsPerDegree, ...
             'flipPSFUpsideDown', flipPSFUpsideDown, ...
             'upsampleFactor',  upsampleFactor, ...
             'noLCA', noLCA, ...
             'name', 'custom human optics');
    
    % Generate the OI from the wavefront map
    theOptics = wvf2oiSpecial(theWVF, micronsPerDegree, pupilDiamMM);
    
end


function [coneLocationsMicrons, LconeIndices, MconeIndices, SconeIndices] = ...
    extractConeLocationsAndTypes(theData)

    % Cone positions
    coneLocations = cell2mat(theData(:,1:2));

    % Cone types
    coneTypes = theData(:,3);
    LconeIndices = find([coneTypes{:}] == 'L');
    MconeIndices = find([coneTypes{:}] == 'M');
    SconeIndices = find([coneTypes{:}] == 'S');

    % Valid cone indices (i.e., not NotClassified)
    validConeIndices = union(LconeIndices, MconeIndices);
    validConeIndices = union(validConeIndices, SconeIndices);

    coneLocationsMicrons = coneLocations(validConeIndices,:);
    coneTypes = coneTypes(validConeIndices);

    % Valid cone locations (non-NaN)
    validConeLocationIndices = find(...
        (~isnan(squeeze(coneLocationsMicrons(:,1)))) & ...
        (~isnan(squeeze(coneLocationsMicrons(:,2)))));

    coneLocationsMicrons = coneLocationsMicrons(validConeLocationIndices,:);
    coneTypes = coneTypes(validConeLocationIndices);

    LconeIndices = find([coneTypes{:}] == 'L');
    MconeIndices = find([coneTypes{:}] == 'M');
    SconeIndices = find([coneTypes{:}] == 'S');
end


function [theConeMosaicData, mosaicMetaData, ZernikeData] = ...
    parseSingleMosaicData(allMosaicsData, iMosaic)

    dataSizeColumns = 7;
    theColumnIndices = (iMosaic-1)*(dataSizeColumns+1) + (1:dataSizeColumns);
    theData = allMosaicsData(:, theColumnIndices);

    [coneLocationsMicrons, LconeIndices, MconeIndices, SconeIndices] = ...
        extractConeLocationsAndTypes(theData);

    % Optics
    ZernikeCoeffNames = squeeze(theData(:,4));
    ZernikeCoeffValues = cell2mat(squeeze(theData(:,5)));

    validZernikeIndices = find(~isnan(ZernikeCoeffValues));
    ZernikeData.coeffNames = ZernikeCoeffNames(validZernikeIndices);
    ZernikeData.coeffValues = ZernikeCoeffValues(validZernikeIndices);

    % Metadata
    paramNames = theData(:,6);
    paramValues = theData(:,7);

    mosaicMetaData = struct(...
        'subjectID', '', ...
        'subjectAge', nan, ...
        'eye', '', ...
        'axialLengthMM', nan, ...
        'horizontalMeridian', '', ...
        'xyEccDegs', nan, ...
        'xyEccMMs', nan, ...
        'retinalMagnificationMicronsPerDegree', nan, ...
        'fovMM', nan, ...
        'LMratio', nan, ...
        'SconePercentage', nan, ...
        'LconeDensityConesPerMM2', nan, ...
        'MconeDensityConesPerMM2', nan, ...
        'SconeDensityConesPerMM2', nan, ...
        'numSelectedCones', nan, ...
        'numNonClassifiedCones', nan, ...
        'LconeString', '', ...
        'MconeString', '', ...
        'SconeString', '', ...
        'nonClassifiedConeString', '', ...
        'coneLocationOrigin', '', ...
        'ZernikeCoeffPupilDiameterMM', nan, ...
        'ZernikeCoeffMeasuredWavelengthNM', nan, ...
        'ZernikeCoeffOptimizeddWavelengthNM', nan ...
     );


     for iLine = 1:numel(paramNames)
        if (find(strcmpi(paramNames{iLine},'subject id')))
            mosaicMetaData.subjectID = paramValues{iLine};
        end
        if (find(strcmpi(paramNames{iLine},'age (years)')))
            mosaicMetaData.subjectAge = paramValues{iLine};
        end

        if (find(strcmpi(paramNames{iLine},'eye')))
            mosaicMetaData.eye = paramValues{iLine};
        end

        if (find(strcmpi(paramNames{iLine},'axial length (mm)')))
            mosaicMetaData.axialLengthMM = paramValues{iLine};
        end

        if (find(strcmpi(paramNames{iLine},'meridian')))
            mosaicMetaData.meridian = paramValues{iLine};
        end

        if (find(strcmpi(paramNames{iLine},'eccentricity (x,y) (deg)')))
            tmp = paramValues{iLine};
            idx = strfind(tmp, ',');
            tmp = cell2mat(tmp);
            mosaicMetaData.xyEccDegs(1) = str2double(tmp(2:idx-1));
            mosaicMetaData.xyEccDegs(2) = str2double(tmp(idx+1:numel(tmp)-1));
        end

        if (find(strcmpi(paramNames{iLine},'eccentricity (x,y) (mm)')))
            tmp = paramValues{iLine};
            idx = strfind(tmp, ',');
            tmp = cell2mat(tmp);
            mosaicMetaData.xyEccMMs(1) = str2double(tmp(2:idx-1));
            mosaicMetaData.xyEccMMs(2) = str2double(tmp(idx+1:numel(tmp)-1));
        end

        if (find(strcmpi(paramNames{iLine},'fov (mm)')))
            tmp = paramValues{iLine};
            idx = strfind(tmp, 'x');
            tmp = cell2mat(tmp);
            mosaicMetaData.fovMM(1) = str2double(tmp(1:idx-1));
            mosaicMetaData.fovMM(2) = str2double(tmp(idx+1:numel(tmp)));
        end


        if (find(strcmpi(paramNames{iLine},'retinal maginification factor (microns/deg)')))
            mosaicMetaData.retinalMagnificationMicronsPerDegree = str2double(cell2mat(paramValues{iLine}));
        end

        if (find(strcmpi(paramNames{iLine}, 'zernike coeffs pupil diameter (mm)')))
            mosaicMetaData.ZernikeCoeffPupilDiameterMM = str2double(cell2mat(paramValues{iLine}));
        end

        if (find(strcmpi(paramNames{iLine}, 'zernike coeffs measured wavelength (nm)')))
            mosaicMetaData.ZernikeCoeffMeasuredWavelengthNM = str2double(cell2mat(paramValues{iLine}));
        end

        if (find(strcmpi(paramNames{iLine}, 'zernike coeffs optimized wavelength (nm)')))
            mosaicMetaData.ZernikeCoeffOptimizeddWavelengthNM = str2double(cell2mat(paramValues{iLine}));
        end

     end % for iLine

    % Zero center the data
    coneLocationsMicrons = bsxfun(@minus, coneLocationsMicrons, median(coneLocationsMicrons,1));
    
    % The Sebasan coordinates are retinal-space referred. 
    % ISETBio are visual space -- referred. So we have to flip the y-coordinate sign
    coneLocationsMicrons(:,2) = -coneLocationsMicrons(:,2);

    % Revert the x-coordinate sign depending on eye and horizontal meridian
    % This will depend on the way the Sabesan lab horizontal position data are encoded
    if (strcmp(mosaicMetaData.eye, 'OD')) && (strcmp(mosaicMetaData.horizontalMeridian, 'Temporal')) || ...
       (strcmp(mosaicMetaData.eye, 'OS')) && (strcmp(mosaicMetaData.horizontalMeridian, 'Nasal'))
       coneLocationsMicrons(:,1) = -coneLocationsMicrons(:,1);
    end

    theConeMosaicData = struct(...
        'coneLocationsMicronsZeroCentered', coneLocationsMicrons, ...
        'LconeIndices', LconeIndices, ...
        'MconeIndices', MconeIndices, ...
        'SconeIndices', SconeIndices, ...
        'mosaicMetaData', mosaicMetaData);
end

function theHumanConeMosaic = generateISETBioConeMosaic(theConeMosaicData, mosaicMetaData)
    
    coneLocationsMicronsZeroCentered = theConeMosaicData.coneLocationsMicronsZeroCentered;
 
    LconeIndices = theConeMosaicData.LconeIndices;
    MconeIndices = theConeMosaicData.MconeIndices;
    SconeIndices = theConeMosaicData.SconeIndices;

    % Revert the theMosaicRect x-coordinate sign depending on eye and horizontal meridian
    % This will depend on the way the Sabesan lab horizontal position data are encoded
    if (strcmp(mosaicMetaData.eye, 'OD')) && (strcmp(mosaicMetaData.horizontalMeridian, 'Temporal')) || ...
       (strcmp(mosaicMetaData.eye, 'OS')) && (strcmp(mosaicMetaData.horizontalMeridian, 'Nasal'))
       theMosaicRect.xDegs = -mosaicMetaData.xyEccDegs(1);
    else
       theMosaicRect.xDegs = mosaicMetaData.xyEccDegs(1);
    end
    theMosaicRect.yDegs = mosaicMetaData.xyEccDegs(2);

    % Mosaic ecc in microns
    coneMosaicEccMicrons = mosaicMetaData.xyEccMMs * 1e3;

    % Generate the default ISETBio mosaic in order to get an estimate
    % of the mean cone aperture, since it is not currently provided for
    % the human mosaics
    theISETbioConeMosaic = cMosaic(...
        'eccentricityDegs', [theMosaicRect.xDegs theMosaicRect.yDegs]);

    meanConeApertureDiameterMicrons = median(theISETbioConeMosaic.coneApertureDiametersMicrons);

    % Add the mosaic ecc to the zero-centered cone locations
    coneLocationsMicrons = bsxfun(@plus, coneLocationsMicronsZeroCentered, coneMosaicEccMicrons);

    conesNum = size(coneLocationsMicrons,1);
    coneApertureDiametersMicrons = ones(1, conesNum)*meanConeApertureDiameterMicrons;

    coneTypes = zeros(1, conesNum);
    assert(conesNum == numel(LconeIndices)+numel(MconeIndices)+numel(SconeIndices), 'L+M+S not equal to total # of cones');
    coneTypes(LconeIndices) = cMosaic.LCONE_ID;
    coneTypes(MconeIndices) = cMosaic.MCONE_ID;
    coneTypes(SconeIndices) = cMosaic.SCONE_ID;

    % Generate struct with custom cone data
    theCustomConeDataStruct = struct(...
        'positionUnits', 'microns', ...
        'positions', coneLocationsMicrons, ...
        'types', coneTypes,...
        'lightGatheringApertureDiameters', coneApertureDiametersMicrons ...                      
        );

    theHumanConeMosaic = cMosaic(...
        'coneData', theCustomConeDataStruct, ...
        'micronsPerDegree', mosaicMetaData.retinalMagnificationMicronsPerDegree);
end


function allSubjectData = importAllSubjectData(filename)

    % All the data
    dataLines = [2, Inf];

    % Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 23);

    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["ConeXLocationmicrons", "ConeYLocationmicrons", "ConeSpectralType", "ZernikeCoeffs", "Values", "Parameter_Name", "Values1", "VarName8", "ConeXLocationmicrons1", "ConeYLocationmicrons1", "ConeSpectralType1", "ZernikeCoeffs1", "Values2", "Parameter_Name1", "Values3", "VarName16", "ConeXLocationmicrons2", "ConeYLocationmicrons2", "ConeSpectralType2", "ZernikeCoeffs2", "Values4", "Parameter_Name2", "Values5"];
    opts.VariableTypes = ["double", "double", "categorical", "string", "double", "string", "string", "string", "double", "double", "categorical", "string", "double", "string", "string", "string", "double", "double", "categorical", "string", "double", "string", "string"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["ZernikeCoeffs", "Parameter_Name", "Values1", "VarName8", "ZernikeCoeffs1", "Parameter_Name1", "Values3", "VarName16", "ZernikeCoeffs2", "Parameter_Name2", "Values5"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["ConeSpectralType", "ZernikeCoeffs", "Parameter_Name", "Values1", "VarName8", "ConeSpectralType1", "ZernikeCoeffs1", "Parameter_Name1", "Values3", "VarName16", "ConeSpectralType2", "ZernikeCoeffs2", "Parameter_Name2", "Values5"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["ConeXLocationmicrons", "ConeYLocationmicrons", "Values", "ConeXLocationmicrons1", "ConeYLocationmicrons1", "Values2", "ConeXLocationmicrons2", "ConeYLocationmicrons2", "Values4"], "ThousandsSeparator", ",");

    % Import the data
    allSubjectData = readtable(filename, opts);

    % Convert to output type
    allSubjectData = table2cell(allSubjectData);
end

