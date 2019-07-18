% COLOR
%
% Color and energy conversion functions
%
% Files
%   cct                           - Calculate correlated color temperature from CIE uv coordinates
%   cct2sun                       - Correlated color temperature to daylight SPD at specified wavelengths
%   chromaticity                  - Compute CIE chromaticity (xy) coordinates from XYZ data
%   chromaticityPlot              - Draw points superimposed on an xy chromaticity diagram
%   colorBlockMatrix              - Create a matrix to render spd data into RGB
%   daylight                      - Generate a daylight SPD with a correlated color temperature (cct)
%   Energy2Quanta                 - Convert energy (watts or joules) to photon units
%   ieLuminanceFromEnergy         - Calculate luminance and related quantities from spectral energy
%   ieLuminanceFromPhotons        - Converts photons into energy and then calls ieLuminanceFromEnergy.
%   ieResponsivityConvert         - Convert sensory responsivity in photons to energy, or vice versa
%   ieScotopicLuminanceFromEnergy - Compute the scotopic (rod) luminance from spectral energy
%   ieSpectraSphere               - Calculate spectra that produce XYZ in a sphere around spectrumE
%   ieXYZFromEnergy               - CIE XYZ values from spectral radiance(w/nm/sr/m2) or irradiance(w/nm/m2)
%   ieXYZFromPhotons              - Convert spectral power distribution in photon units into CIE XYZ
%   initDefaultSpectrum           - Create a wavelength spectrum structure and attach it to an ISET object
%   Quanta2Energy                 - Convert quanta (photons) to energy (watts)
%   RGB2XWFormat                  - Transform an RGB form matrix into an XW (space-wavelength) matrix
%   spd2cct                       - Convert a spectral power distribution to a correlated color temperature
%   XW2RGBFormat                  - Convert XW format data to RGB format
