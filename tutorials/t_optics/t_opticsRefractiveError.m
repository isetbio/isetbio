% Examine the effect of refractive error the on retinal image
%
% Description:
%    Generates an ISETBio scene from a jpg file. Passes the scene via
%    optics in which the user controls the additional refractive error.
%

% History:
%    08/24/21       NPC Wrote it. Copyright, ISETBIO Team, 2021
%    01/12/23       NPC Added an extensive comparison

function t_opticsRefractiveError

    % Your favorite image source
    fname = 'ma_blkbkjackal_412.jpg';
    
    % Mean luminance of the scene
    meanLuminanceCdPerM2 = 100;
    
    % Field of view of the scene
    fieldOfViewDegs = 2;
    
    % Generate an ISETBIo hyperspectral scene assuming it is displayed on a
    % typical Apple LCD display
    
    scene = sceneFromFile(fname, 'rgb', meanLuminanceCdPerM2, 'LCD-Apple.mat');
    scene = sceneSet(scene, 'wAngular', fieldOfViewDegs);
    
    % Optics at a given eccentricity
    eccentricityDegs = [0 0];
    
    % Pupil diameter in MM
    pupilDiameterMM = 4.0;
    
    % Human retina magnification factor (retinal microns/deg of visual angle)
    retinalMagnificationMicronsPerDegree = 300;
    
    % How much refractive error to include
    additionalRefractiveErrorDiopters = 0.75;
    
    % Compare effect of adding an additional refractive error for one Artal
    % subject (#11). Also get the Z coeffs of the optics without and with
    % the additional refractive errror and compate them
    ArtalSubjectID = 11;
    [theZernikeCoeffsWithoutAdditionalRefractiveError, ...
     theZernikeCoeffsWithAdditionalRefractiveError] = compareEffectOfAdditionalRefractiveError(...
        'Artal', ArtalSubjectID, additionalRefractiveErrorDiopters, ...
        eccentricityDegs, pupilDiameterMM, retinalMagnificationMicronsPerDegree, ...
        scene);

    fprintf('The Zcoeffs of the generated optics (with/without additional refactive error) differ in the following indices: %d\n', ...
       find(theZernikeCoeffsWithoutAdditionalRefractiveError ~= theZernikeCoeffsWithAdditionalRefractiveError))
    fprintf('The defocus Zernike coefficient is Z(%d)\n', wvfOSAIndexToVectorIndex('defocus'));



    % Compare effect for one Polans subject (#2)
    PolansSubjectID = 2;
    [theZernikeCoeffsWithoutAdditionalRefractiveError, ...
     theZernikeCoeffsWithAdditionalRefractiveError] = compareEffectOfAdditionalRefractiveError(...
        'Polans', PolansSubjectID, additionalRefractiveErrorDiopters, ...
        eccentricityDegs, pupilDiameterMM, retinalMagnificationMicronsPerDegree, ...
        scene);


    fprintf('The Zcoeffs of the generated optics (with/without additional refactive error) differ in the following indices: %d\n', ...
       find(theZernikeCoeffsWithoutAdditionalRefractiveError ~= theZernikeCoeffsWithAdditionalRefractiveError))
    fprintf('The defocus Zernike coefficient is Z(%d)\n', wvfOSAIndexToVectorIndex('defocus'));



    % Compare effect for diffraction-limited optics
    PolansSubjectID = 0;
    [theZernikeCoeffsWithoutAdditionalRefractiveError, ...
     theZernikeCoeffsWithAdditionalRefractiveError] = compareEffectOfAdditionalRefractiveError(...
        'Polans', PolansSubjectID, additionalRefractiveErrorDiopters, ...
        eccentricityDegs, pupilDiameterMM, retinalMagnificationMicronsPerDegree, ...
        scene);


     fprintf('The Zcoeffs of the generated optics (with/without additional refactive error) differ in the following indices: %d\n', ...
       find(theZernikeCoeffsWithoutAdditionalRefractiveError ~= theZernikeCoeffsWithAdditionalRefractiveError))
    fprintf('The defocus Zernike coefficient is Z(%d)\n', wvfOSAIndexToVectorIndex('defocus'));
    
end


function [theZernikeCoeffsWithoutAdditionalRefractiveError, theZernikeCoeffsWithAdditionalRefractiveError] = ...
    compareEffectOfAdditionalRefractiveError(opticsDataBase, opticsSubjectID, additionalRefractiveErrorDiopters, ...
    eccentricityDegs, pupilDiameterMM, retinalMagnificationMicronsPerDegree, scene)


    switch (opticsDataBase)
        case 'Artal'
            % Generate optics based on human wavefront aberration measurements and the
            % subject's normal refractive error
            [theDefaultOI, theDefaultRefractiveErrorPSF, ...
             thePSFspatialSupportXminutes, ...
             thePSFspatialSupportYminutes,...
             thePSFwavelengthSupportNanometers, ...
             theZernikeCoeffsWithoutAdditionalRefractiveError] = ArtalOptics.oiForSubjectAtEccentricity(opticsSubjectID, ...
                  'right eye', eccentricityDegs, pupilDiameterMM, ...
                   sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
                  'subtractCentralRefraction', false, ...
                  'refractiveErrorDiopters', 0.0, ...
                  'zeroCenterPSF', true);
        
            % Generate optics based on human wavefront aberration measurements and
            % the passed  additional refractive error
            [theAdditionalRefractiveErrorOI, theAdditionalRefractiveErrorPSF, ...
             thePSFspatialSupportXminutes, ...
             thePSFspatialSupportYminutes,...
             thePSFwavelengthSupportNanometers, ...
             theZernikeCoeffsWithAdditionalRefractiveError] = ArtalOptics.oiForSubjectAtEccentricity(opticsSubjectID, ...
                  'right eye', eccentricityDegs, pupilDiameterMM, ...
                  sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
                  'refractiveErrorDiopters', additionalRefractiveErrorDiopters, ...
                  'zeroCenterPSF', true);

        case 'Polans'
            % Generate optics based on human wavefront aberration measurements and the
            % subject's normal refractive error
            [theDefaultOI, theDefaultRefractiveErrorPSF, ...
             thePSFspatialSupportXminutes, ...
             thePSFspatialSupportYminutes,...
             thePSFwavelengthSupportNanometers, ...
             theZernikeCoeffsWithoutAdditionalRefractiveError] = PolansOptics.oiForSubjectAtEccentricity(opticsSubjectID, ...
                  'right eye', eccentricityDegs, pupilDiameterMM, ...
                   sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
                  'subtractCentralRefraction', false, ...
                  'refractiveErrorDiopters', 0.0, ...
                  'zeroCenterPSF', true);
        
            % Generate optics based on human wavefront aberration measurements and
            % the passed  additional refractive error
            [theAdditionalRefractiveErrorOI, theAdditionalRefractiveErrorPSF, ...
             thePSFspatialSupportXminutes, ...
             thePSFspatialSupportYminutes,...
             thePSFwavelengthSupportNanometers, ...
             theZernikeCoeffsWithAdditionalRefractiveError] = PolansOptics.oiForSubjectAtEccentricity(opticsSubjectID, ...
                  'right eye', eccentricityDegs, pupilDiameterMM, ...
                  sceneGet(scene, 'wave'), retinalMagnificationMicronsPerDegree, ...
                  'refractiveErrorDiopters', additionalRefractiveErrorDiopters, ...
                  'zeroCenterPSF', true);

        otherwise
            error('Unknown optics database: ''%s''.', opticsDataBase)
    end

    % Compute the retinal images for the two OIs
    theDefaultOI = oiCompute(scene,theDefaultOI);
    theAdditionalRefractiveErrorOI = oiCompute(scene,theAdditionalRefractiveErrorOI);

    % Visualize scene and its retinal image under the employed optics
    hFig = figure();
    set(hFig, 'Position', [750 50 1240 1240], 'Name', sprintf(''));
    subplot(3,2,[1 2]);
    image(sceneGet(scene, 'rgb image'));
    axis 'image'
    title('test scene');
    
    subplot(3,2,3);
    image(oiGet(theDefaultOI, 'rgb image'));
    axis 'image'
    title(sprintf('retinal image\ndefault refractive error'));

    subplot(3,2,4);
    % Display PSF at 550 nm
    [~,idx] = min(abs(thePSFwavelengthSupportNanometers-550));
    theDefaultRefractiveErrorPSF = squeeze(theDefaultRefractiveErrorPSF(:,:,idx));
    imagesc(thePSFspatialSupportXminutes, thePSFspatialSupportYminutes,theDefaultRefractiveErrorPSF);
    axis 'image'
    set(gca, 'XLim', [-5 5], 'YLim', [-5 5])
    title(sprintf('PSF\ndefault refractive error'));

    subplot(3,2,5);
    image(oiGet(theAdditionalRefractiveErrorOI , 'rgb image'));
    axis 'image'
    title(sprintf('retinal image\nadditional refractive error: %2.2f D', additionalRefractiveErrorDiopters));
    
    subplot(3,2,6);
    theDefaultRefractiveErrorPSF = squeeze(theAdditionalRefractiveErrorPSF(:,:,idx));
    imagesc(thePSFspatialSupportXminutes, thePSFspatialSupportYminutes,theDefaultRefractiveErrorPSF);
    axis 'image'
    set(gca, 'XLim', [-5 5], 'YLim', [-5 5])
    title(sprintf('PSF\nadditional refractive error'));
end

