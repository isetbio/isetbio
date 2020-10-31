function [sceneLRGBimage, sceneSRGBimage, sceneLMScontrastsImage, sceneLMSexcitationsImage] = ...
    sceneRepresentations(theScene, presentationDisplay)

    emittedRadianceImage = sceneGet(theScene, 'energy');
    displaySPDs = displayGet(presentationDisplay, 'spd');
    displayWavelengths = displayGet(presentationDisplay, 'wave');
    [sceneLRGBimage, sceneSRGBimage] = displayRadianceToDisplayRGB(emittedRadianceImage, displaySPDs);
        
    % Load the 2-deg Stockman cone fundamentals on a wavelength support matching the display
    coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);

    % Compute the LMS cone contrasts of the emitted radiance image
    [sceneLMScontrastsImage, sceneLMSexcitationsImage] = displayRadianceToLMS(emittedRadianceImage, coneFundamentals);
end
