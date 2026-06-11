function visualizeSpectralRadianceSlice(spatialSupport, wavelengthSupport, photonRate, targetPos, sceneName)

% Extract X,Y spatial supports
spatialSupportX = squeeze(spatialSupport(:,:,1));
spatialSupportY = squeeze(spatialSupport(:,:,2));

% Find target pixel
[~,targetPixelIndex] = ...
    min(sum(([spatialSupportX(:) spatialSupportY(:)]-reshape(targetPos, [1 2])).^2,2));
[targetPixelRow, targetPixelCol] = ...
    ind2sub(size(spatialSupportX), targetPixelIndex);

% Extract spectral slice of photons at target pixel
% photonRate is in photons/sec/nm/sr/m^2
photonRateSpectralSlice = squeeze(photonRate(targetPixelRow, targetPixelCol,:));

% Compute energy in Watts/sr/nm/m2 from the photon rate
wavelengBandwidth = wavelengthSupport(2)-wavelengthSupport(1);
energyInWattsPerSrPerNMPerM2 = Quanta2Energy(wavelengthSupport, photonRateSpectralSlice);

% Read-in the human luminosity function, Vlambda
fName = fullfile(isetbioDataPath, 'human', 'luminosity.mat');
Vlambda = ieReadSpectra(fName, wavelengthSupport);

% Read-in the CIE XYZ color matching functions
fName = fullfile(isetbioDataPath, 'human', 'XYZ.mat');
XYZ = ieReadSpectra(fName, wavelengthSupport);

% Compute target luminance by integrating the dot product 
% of energy with VLambda over wavelength
radiometricToPhotometricConversionFactor = 683; % Cd*sr*/watt
targetLuminanceCdPerM2 = radiometricToPhotometricConversionFactor * ...
                  dot(energyInWattsPerSrPerNMPerM2(:),Vlambda(:)) * ...
                  wavelengBandwidth;

% Compute target chromaticity by integrating the dot product 
% of energy with X,Y,Z colormatching functions over wavelength
targetX = dot(energyInWattsPerSrPerNMPerM2(:),XYZ(:,1),1);
targetY = dot(energyInWattsPerSrPerNMPerM2(:),XYZ(:,2),1);
targetZ = dot(energyInWattsPerSrPerNMPerM2(:),XYZ(:,3),1);
targetChroma = [targetX targetY]/(targetX+targetY+targetZ);

% Plot spectral radiance slice
figure(); clf;
plot(wavelengthSupport, photonRateSpectralSlice, 'bo-', ...
    'MarkerFaceColor', [0 0.8 0.9], 'MarkerSize', 12, ...
    'LineWidth', 1.5);
set(gca, 'YLim', [min(photonRate(:)) max(photonRate(:))], ...
    'XLim', [wavelengthSupport(1) wavelengthSupport(end)]);
set(gca, 'FontSize', 16);
xlabel('\it wavelength (nm)');
ylabel('photon rate (photons/sec/nm/sr/m^2)');
title(sprintf('%s\nluminance at (%2.2f %2.2f): %2.2f cd/m2\n chromaticity: (%2.3f %2.3f)', ...
    sceneName, spatialSupportX(targetPixelIndex), spatialSupportY(targetPixelIndex), ...
    targetLuminanceCdPerM2, targetChroma(1), targetChroma(2)));
drawnow;
end
