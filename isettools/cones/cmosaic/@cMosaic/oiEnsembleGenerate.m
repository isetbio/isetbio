function [oiEnsemble, psfEnsemble, zCoeffs, oiSamplingGridDegs, theMergingWeights] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin)
% Create an ensemble of optical images and psfs
%
% NOTE:  BW and DHB/NC should examine the organization of this method and
% simplify it.  
% 
% As written, the method has its own special wavefront calculation tools
% (makeWVF, computePSFandOTF) rather than relying on the core wavefront
% tools (wvfCreate, wvfCompute(), wvfGet(). Consequently, the number of
% parameters and complexity of the methods is much greater than we would
% like (IMHO).  Also, there are many special numbers in here for specifying
% degrees, and the particulars of this call.  That means as we improve the
% general functions, and validate and test them, we do not similarly
% improve this function.  An example is the new addition of a pupil
% aperture function. This format does not have such a function, which might
% be useful for modeling scratched corneas.  Or oddly shaped pupils.  Or
% ...
%
% I think we want the call to the subject in any database to simply return
% the Zernike coefficients in a wavefront structure.  Then we make a few
% calls to wvfCompute(), wvfGet(this and that), and we are done.
%
%
% Brief description
%  The varargin parameters specify the dataset and other processing
%  parameters.
%
% Input:
%   obj:  A cMosaic
%
% Optional key/val pairs
%  zernikeDataBase           - Either one of {'Polans2015', 'Artal2012'}
%  subjectID                 - Subject number in the data base. Default 6.
%  pupilDiameterMM           - Pupil diameter in millimeters (default 3mm)
%  wavefrontSpatialSamples   - Default is 301 samples
%  subtractCentralRefraction - Fixes up the data a bit. Default false.
%  zeroCenterPSF             - Default true
%  deNoisedZernikeCoefficients - Seems deprecated to me (BW).
%  flipPSFUpsideDown         - Default true (aligns with image)
%  defocusMicrons            - Add this value to defocus zcoeff
%
% Outputs:
%   oiEnsemble          -  Cell array of OIs
%   psfEnsemble         -  Cell array of PSFs
%   zCoeffs             -  Zernike polynomial coeffcieints
%   oiSamplingGridDegs  -  the OI sampling grid in degs
%
% Description:
%   In which we explain more about the processing parameters.  NC and BW to
%   do together.
%
% See also
%   cMosaic (main class)
%

% History
%   11/22/22  dhb  Add Marimont/Wandell optics options as best I could.

% Help
if (ischar(oiSamplingGridDegs)) && (strcmp(oiSamplingGridDegs,'help'))
    doc('cMosaic.oiEnsembleGenerate');
    return;
end

% Parse input
p = inputParser;
p.addRequired('obj', @(x)(isa(x,'cMosaic')));
p.addRequired('oiSamplingGridDegs',  @(x)( (isstruct(x)) || (isnumeric(x) && (size(x,2) == 2) && (all(isreal(x)))) ) );
p.addParameter('zernikeDataBase', 'Polans2015', @(x)(ismember(x, {'Polans2015', 'Artal2012', 'MarimontWandell', 'Thibos2002'})));
p.addParameter('warningInsteadOfErrorForBadZernikeCoeffs', false, @islogical);
p.addParameter('subjectID', 6, @isscalar);
p.addParameter('pupilDiameterMM', 3.0, @isscalar);
p.addParameter('inFocusWavelength', 550, @isscalar);
p.addParameter('wavefrontSpatialSamples', 301, @isscalar);
p.addParameter('subtractCentralRefraction', false, @islogical);
p.addParameter('zeroCenterPSF', true, @islogical);
p.addParameter('flipPSFUpsideDown', true, @islogical);
p.addParameter('visualizedSamplingGrid', false, @islogical);
p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
p.addParameter('noLCA',false,@islogical);
p.parse(obj, oiSamplingGridDegs, varargin{:});

oiSamplingGridDegs = p.Results.oiSamplingGridDegs;
zernikeDataBase = p.Results.zernikeDataBase;
pupilDiamMM = p.Results.pupilDiameterMM;
inFocusWavelength = p.Results.inFocusWavelength;
subjectID = p.Results.subjectID;
subtractCentralRefraction = p.Results.subtractCentralRefraction;
wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
zeroCenterPSF = p.Results.zeroCenterPSF;
flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
upSampleFactor = p.Results.upsampleFactor;
warningInsteadOfErrorForBadZernikeCoeffs = p.Results.warningInsteadOfErrorForBadZernikeCoeffs;
visualizedSamplingGrid = p.Results.visualizedSamplingGrid;

if (isstruct(oiSamplingGridDegs))
    theOIsampligGridStruct = oiSamplingGridDegs;

    % Generate sampling grid from passed struct
    oiSamplingGridDegs = oiSamplingGridFromStruct(theOIsampligGridStruct);

    % Compute weight with which we will merge the cone mosaic responses to
    % the different OIs
    theMergingWeights = computeMergingWeights(obj, ...
        oiSamplingGridDegs, theOIsampligGridStruct.spacingDegs, theOIsampligGridStruct.weightingType);
else
    theMergingWeights = [];
end




% Generate the oiEnsemble
oiNum = size(oiSamplingGridDegs,1);
oiEnsemble = cell(1, oiNum);
psfEnsemble = cell(1, oiNum);

switch (zernikeDataBase)
    case 'MarimontWandell'
        % This doesn't use wavefront optics but implements the Marimont and
        % Wandell optics model.  This is for foveal viewing but at some
        % risk we allow it for any eccentricity but warn the user with
        % a printout.

        for oiIndex = 1:oiNum
            %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
            %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);

            if (targetEcc(1) ~= 0 || targetEcc(2) ~= 0)
                fprintf(2,'Marimont/Wandell optics not available off the fovea. Computing for hEcc = 0 and vEcc = 0\n');
                targetEcc(1) = 0;
                targetEcc(2) = 0;
            end

            % Create Marimont/Wandell oi.  I don't see how to set
            % microns per degree easily here, and so we ignore the microns
            % per degree property of the passed mosaic object.
            [theOI] = oiCreate('human',pupilDiamMM,'wave',obj.wave);
            theOptics = oiGet(theOI,'optics');
            thePSF = opticsGet(theOptics,'psf data');

            % God save us, the comments in oiGet don't provide the unit
            % options and say that we aren't sure whether this comes back
            % as X/Y or Y/X. Looking through the code, I think passing
            % units of frequency to be cycles per degree and thus the units
            % of psf support to be degrees.
            psfSupport = opticsGet(theOptics,'psf support','cyclesperdeg');
            psfSupportMinutesX = psfSupport{1}*60;
            psfSupportMinutesY = psfSupport{2}*60;
            psfSupportWavelength = opticsGet(theOptics,'wave');
            zCoeffs = [];

            if (isempty(theOI))
                if (warningInsteadOfErrorForBadZernikeCoeffs)
                    fprintf(2,'Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye\n', obj.whichEye, subjectID);

                    oiEnsemble = [];
                    psfEnsemble = [];
                    zCoeffs = [];
                    return;
                else
                    error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
                end

            end

            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', []);
        end

    case 'Artal2012'
        % Looks like Artal optics now accepts refractive error in diopters.
        % Commented out this warning. DHB.
        %
        % % Make sure refractive error is zero, because Artal version of
        % % oiForSubjectAtEccentricity doesn't understand the
        % % 'refractiveErrorMicrons' key/value pair.
        % if (p.Results.refractiveErrorDiopters ~= 0)
        %     error('Artal optics does not currently accept refractiveErrorDiopters key/value pair');
        % end

        % Artal optics
        for oiIndex = 1:oiNum
            %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
            %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);

            if (abs(targetEcc(2)) > 0.1)
                fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
                targetEcc(2) = 0;
            end

            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                obj.whichEye, targetEcc, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                'inFocusWavelength', inFocusWavelength, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'upsampleFactor', upSampleFactor, ...
                'noLCA',p.Results.noLCA, ...
                'refractiveErrorDiopters', p.Results.refractiveErrorDiopters);

            if (isempty(theOI))
                if (warningInsteadOfErrorForBadZernikeCoeffs)
                    fprintf(2,'Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye\n', obj.whichEye, subjectID);

                    oiEnsemble = [];
                    psfEnsemble = [];
                    zCoeffs = [];
                    return;
                else
                    error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
                end

            end

            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', zCoeffs);
        end

    case 'Polans2015'
        % Polans optics
        for oiIndex = 1:oiNum
            fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
                zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);

            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
                PolansOptics.oiForSubjectAtEccentricity(subjectID, ...
                obj.whichEye, targetEcc, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                'inFocusWavelength', inFocusWavelength, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'upsampleFactor', upSampleFactor, ...
                'noLCA',p.Results.noLCA, ...
                'refractiveErrorDiopters', p.Results.refractiveErrorDiopters);

            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', zCoeffs);
        end

    case 'Thibos2002'

        idx = find(ThibosOptics.constants.availableMeasurementPupilDiamsMM >= pupilDiamMM);
        if (isempty(idx))
            error('pupilDiameterMM must be <= %f', max(ThibosOptics.constants.availableMeasurementPupilDiamsMM));
        else
            measurementPupilDiameterMM =  ThibosOptics.constants.availableMeasurementPupilDiamsMM(idx(1));
        end

        % Thibos optics
        for oiIndex = 1:oiNum
            %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
            %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);
            if (abs(targetEcc(1)) > 0.1 || abs(targetEcc(2)) > 0.1)
                fprintf(2,'Thibos optics not available off the fovea. Computing for hEcc = 0 and vEcc = 0\n');
            end

            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
                ThibosOptics.oiForSubjectAtEccentricity(subjectID, ...
                obj.whichEye, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                'inFocusWavelength', inFocusWavelength, ...
                'measurementPupilDiameterMM', measurementPupilDiameterMM, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'upsampleFactor', upSampleFactor, ...
                'noLCA',p.Results.noLCA, ...
                'refractiveErrorDiopters', p.Results.refractiveErrorDiopters);

            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', zCoeffs);
        end

    otherwise
        error('Unknown Zernike database specified');

end

if (~isempty(theMergingWeights)) && (visualizedSamplingGrid)

    % Visualize the merging weights
    for oiPosIndex = 1:size(theMergingWeights,1)
        visualizeTheMergingWeights(obj, theMergingWeights(oiPosIndex,:), oiSamplingGridDegs);
    end
    visualizeThePSFs(obj, psfEnsemble, oiSamplingGridDegs)
end

end


function oiSamplingGridDegs = oiSamplingGridFromStruct(oiSamplingGridStruct)
    % Generate a regular hexagonal grid with specified 
    % spacing, height, width at the specified eccentricity 

    spacing = oiSamplingGridStruct.spacingDegs;
    width = oiSamplingGridStruct.widthDegs;
    height = oiSamplingGridStruct.heightDegs;
    centerCoords = oiSamplingGridStruct.eccentricityDegs;

    dx = spacing;
    dy = spacing * sqrt(3) / 2;
    
    % Number of columns and rows
    colsNum = ceil(width  / dx) + 2;
    rowsNum = ceil(height / dy) + 2;
    
    x = []; y = [];
    for row = 0:(rowsNum - 1)
        x_offset = mod(row, 2) * (dx / 2);
        for col = 0 : colsNum - 1
            xi = col * dx + x_offset;
            yi = row * dy;
            if (xi <= width) && (yi <= height)
                x(end+1) = xi;  
                y(end+1) = yi;  
            end
        end
    end

    x = x - mean(x(:));
    y = y - mean(y(:));
    oiSamplingGridDegs = bsxfun(@plus, [x(:) y(:)], centerCoords);

    % Sort them according to their ecc
    d = sum(oiSamplingGridDegs.^2,2);
    [~,idx] = sort(d, 'ascend');
    oiSamplingGridDegs = oiSamplingGridDegs(idx,:);
end

function mergingWeights = computeMergingWeights(theConeMosaic, oiSamplingGridDegs, spacingDegs, weightingType)

    conesNum = theConeMosaic.conesNum;
    opticsSamplingPositionsNum = size(oiSamplingGridDegs,1);
    theConeOpticsSamplingPositionWeights = zeros(opticsSamplingPositionsNum, conesNum);

    minDistanceDegs = 0.01;

    if (strcmp(weightingType, 'Gaussian'))
        sigmaDegs = 0.5*spacingDegs;
    end

    % Compute weights of all cones to all oiPositions
    for oiPos = 1:opticsSamplingPositionsNum
        distanceDegs = sqrt(sum((bsxfun(@minus,  theConeMosaic.coneRFpositionsDegs,oiSamplingGridDegs(oiPos,:))).^2,2));
        % Avoid dividing with too small distances
        distanceDegs = max(minDistanceDegs, distanceDegs);

        if (strcmp(weightingType, 'Gaussian'))
            theWeights = exp(-0.5*(distanceDegs/sigmaDegs).^2);
        elseif (strcmp(weightingType, 'raised cosine'))
            theWeights = 0.5 * (1 + cos(pi * min(1,distanceDegs/spacingDegs)));
        else
            error('Unknown weightingType: ''%s''.', weightingType);
        end

        theConeOpticsSamplingPositionWeights(oiPos,:) = theWeights;
    end

    % Sum of all weights of a single cone to all oiPositions
    s = sum(theConeOpticsSamplingPositionWeights,1);

    % Ensure that the sum of all weights of a single cone to all
    % oiPositions is equal to 1.0;
    mergingWeights = bsxfun(@times, theConeOpticsSamplingPositionWeights, 1./s);
    s = sum(mergingWeights,1);
    assert((abs(min(s)-1) < 10*eps) || (abs(max(s)-1) < 10*eps), 'error: weights do not sum to 1');
end



function  hFig = visualizeThePSFs(theConeMosaic, psfEnsemble, oiSamplingGridDegs)

    opticsSamplingPositionsNum = size(oiSamplingGridDegs, 1);

    hFig = figure(3002); clf;
    if (theConeMosaic.sizeDegs(1)>theConeMosaic.sizeDegs(2))
        % wide mosaic
        figureWidthPixels = 1200;
        aspectRatio = theConeMosaic.sizeDegs(1)/theConeMosaic.sizeDegs(2);
        figureHeightPixels = round(figureWidthPixels/aspectRatio);
    else
        figureHeightPixels = 1200;
        aspectRatio = theConeMosaic.sizeDegs(2)/theConeMosaic.sizeDegs(1);
        figureWidthPixels = round(figureHeightPixels/aspectRatio);
    end

    set(hFig, 'Position', [10 10 figureWidthPixels figureHeightPixels], 'Color', [1 1 1]);

    for oiPos = 1:opticsSamplingPositionsNum
        
        width = 0.2;
        axPosition(1) = 0.5*(1+(oiSamplingGridDegs(oiPos,1) - theConeMosaic.eccentricityDegs(1))/(0.52*theConeMosaic.sizeDegs(1)))-0.5*width;
        axPosition(2) = 0.5*(1+(oiSamplingGridDegs(oiPos,2) - theConeMosaic.eccentricityDegs(2))/(0.52*theConeMosaic.sizeDegs(2))) - 0.5*width;
        axPosition(3:4) = width;

        ax = axes('Position', axPosition);
        set(ax, 'Color', [0 0 0]);
        thePSF = psfEnsemble{oiPos};
        [~, wIdx] = min(abs(thePSF.supportWavelength-550));
        wavePSF = squeeze(thePSF.data(:,:,wIdx));
        zLevels = 0.1:0.1:0.9;
        % half a degree
        xyRangeArcMin = 15*[-1 1];
        PolansOptics.renderPSF(ax, ...
            thePSF.supportX, thePSF.supportY, wavePSF/max(wavePSF(:)), ...
            xyRangeArcMin, zLevels,  gray(1024), [0 0 0], ...
            'plotTitle',  sprintf('%2.1f,%2.1f', oiSamplingGridDegs(oiPos,1), oiSamplingGridDegs(oiPos,2)));
        box(ax, 'on')
        set(ax, 'XColor', 0.2*[1 1 1], 'YColor', 0.2*[1 1 1]);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        xlabel(ax, '');
        ylabel(ax, '');
        colormap(ax, gray(1024));
        drawnow;
    end
    pause(0.5)
end

function hFig = visualizeTheMergingWeights(theConeMosaic, theMergingWeights, oiSamplingGridDegs)
    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 700 1300], 'Color', [1 1 1]);
    ax = subplot('Position', [0.05 0.05 0.94 0.94]);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', reshape(theMergingWeights, [1 1 numel(theMergingWeights)]), ...
        'activationRange', [0 1], ...
        'plotTitle', ' ');
    hold(ax, 'on');
    plot(ax, oiSamplingGridDegs(:,1), oiSamplingGridDegs(:,2), 'wx', 'MarkerSize', 12, 'LineWidth', 1.5);
end
