function fName = ieSaveSIDataFile(psf, wave, umPerSamp, fName)
% Write file with data for shift-invariant optics
%
% Syntax:
%   fName = ieSaveSIDataFile(psf, wave, umPerSamp, [fName])
%
% Description:
%    The shift-invariant optics uses one point spread for every wavelength.
%    We save the data used to create optics in this file. We then create
%    the optics structure from these data using siSynthetic.
%
% Inputs:
%    psf       - row x col x wavelength array containing the psfs
%    wave      - list of wavelengths
%    umPerSamp - spacing (in microns) between samples in the psf data
%    fName     - (Optional) output file name. Selected when not provided.
%
% Outputs:
%    fName     - The filename the SI Data File has been saved under.
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010
%    11/27/17  jnm  Formatting
%    11/29/17  jnm  Added note in notes section, and note in example
%

% Examples:
%{
    % [Note: JNM - Opens GUI to retrieve file as fName is not specified.]
    psf = rand(128, 128, 31);
    wave = 400:10:700;
    umPerSamp = [0.25, 0.25];
    fName = ieSaveSIDataFile(psf, wave, umPerSamp)
 
    oi = vcGetObject('OI');
    optics = siSynthetic('custom', oi, fName);  
%}

if notDefined('psf'), error('psf volume required'); end
if notDefined('wave'), error('wavelength samples required (nm)'); end
if notDefined('umPerSamp'), error('Microns per sample required'); end
if notDefined('fName'), fName = vcSelectDataFile('stayPut', 'w'); end

notes.timeStamp = datestr(now);

save(fName, 'psf', 'wave', 'umPerSamp', 'notes');

end