function [optics, inName, outName] = ...
    siConvertRTdata(inName, fieldHeight, outName)
% Convert RT optics to a custom optics structure data from one field height
%
% Syntax:
%   optics = siConvertRTdata([inName], [fieldHeight], [outName])
%
% Description:
%    fieldHeight specified in meters
%
%    The RT data are a good source of examples for single PSFs from real
%    lenses.  This routine creates an optic structure used for
%    shift-invariant calculations but with an OTF/PSF drawn from one of the
%    field heights in a ray trace (Zemax) calculation.
%
%    The saved file can then be used for shift-invariant calculations with
%    a custom calculation.
%
% Inputs:
%    inName      - (Optional) String. The input filename. Default queries
%                  user to select via GUI.
%    fieldHeight - (Optional) Numeric. The field height. Default queries
%                  user to input via GUI. Query is for mm, but variable is
%                  stored in meters.
%    outName     - (Optional) String. The output file name. Default queries
%                  user to select via GUI.
%
% Outputs:
%    optics      - Struct. The optics structure.
%    inName      - String. The input filename.
%    outName     - String. The output filename.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/21/17  dhb  Replace fftshift(fft2(psf) with call to PsfToOtf with
%                   an ifftshift.  This is how I think it should be done.
%                   Because this routine is not currently called anywhere
%                   in isetbio, I don't have a good way to test. But the
%                   change matches the change I implemented and tested in
%                   siSynthetic.
%    03/16/18  jnm  Formatting.
%    04/07/18  dhb  Comment wrt broken example and what to do about it.

% Examples:
%{
    % ETTBSkip.  Example is broken.  And it requires user input.
    % Should try to fix so it runs, then leave in the ETTBSkip because
    % it requires user input. I think the problem is files it wants
    % to read are not where the example expects them, and that it isn't
    % clear in the case where it prompts for input what file we need to
    % navigate to and select.

    baseDir = fullfile(isetbioDataPath, 'optics');
    inName = fullfile(baseDir, 'rtZemaxExample.mat');

    siConvertRTdata;

    fieldHeight = 0.5;
    siConvertRTdata(inName, fieldHeight, ...
        fullfile(baseDir, 'siZemaxExample05.mat'));

    fieldHeight = 1.0;
    siConvertRTdata(inName, fieldHeight, ...
        fullfile(baseDir, 'siZemaxExample10.mat'));
%}

if notDefined('inName'), inName = vcSelectDataFile; end
if notDefined('fieldHeight')
    fieldHeight = ieReadNumber('Enter field height (mm)', 0, '%.02f');
    fieldHeight = fieldHeight / 1000;  % fieldHeight must be in meters
end

% Read in the ray traced optics file
tmp = load(inName);
rtOptics = tmp.optics;
clear tmp;

% Figure out the nyQuist
rtWave = opticsGet(rtOptics, 'rtPSFWavelength');
dx = opticsGet(rtOptics, 'rtPSFSpacing', 'm');
rtSupport = opticsGet(rtOptics, 'rtSupport', 'm');
nSamples = size(rtSupport, 1);

nyquistF = 1 ./ (2 * dx);  % Line pairs (cycles) per meter

OTF = zeros(nSamples, nSamples, length(rtWave));
for ii = 1:length(rtWave)
    psf = opticsGet(rtOptics, 'rtpsfdata', fieldHeight, rtWave(ii));
    psf = psf / sum(psf(:));

    % Use PsfToOtf to make the change, and then put center in upper right
    % to match isetbio conventions.  Commented out below is the older code,
    % which may or may not do the same thing
    [~, ~, centeredOTF] = PsfToOtf([], [], psf);
    OTF(:, :, jj) = ifftshift(centeredOTF);
    % OTF(:, :, ii) = fftshift(fft2(psf));
end

% Check this - we converted from mm to meters ... make sure everything
% plots and looks OK
fx = unitFrequencyList(nSamples) * nyquistF(2);
fy = unitFrequencyList(nSamples) * nyquistF(1);
% [FY, FX] = meshgrid(fy, fx);
% figure;
% mesh(FY, FX, abs(OTF(:, :, ii)))
% figure;
% mesh(FY, FX, OTF(:, :, ii))

optics = opticsCreate;

% We may have a problem with the meters scale here ... we were probably
% setting cyc/mm and now we are setting cyc/meter ...
% April 26, 2008
optics = opticsSet(optics, 'otffunction', 'custom');
optics = opticsSet(optics, 'otfData', OTF);
optics = opticsSet(optics, 'otffx', fx);
optics = opticsSet(optics, 'otffy', fy);
optics = opticsSet(optics, 'otfwave', rtWave);

if notDefined('outName'), outName = vcSelectDataFile('stayput', 'w'); end
vcSaveObject(optics, outName);

end