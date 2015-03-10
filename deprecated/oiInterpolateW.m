function oi = oiInterpolateW(oi,newWave)
%Wavelength interpolation for optical image data
%
%  oi = oiInterpolateW(oi,[newWave])
%
% Purpose:
%    Interpolate the wavelength dimension of an optical image.
%   
% Examples:
%   oi = oiInterpolateW(oi,[400:10:700])
%
% Copyright ImagEval Consultants, LLC, 2003.

% Programming ... usually we interpolate in the scene, not the oi.
disp('Warning.  oiInterpolateW not debugged yet.')

if notDefined('oi'), [val,oi] = vcGetSelectedObject('oi'); end
handles = ieSessionGet('opticalimagehandle');

% Note the current oi properties
row = oiGet(oi,'row'); 
col = oiGet(oi,'col'); 
nWave = oiGet(oi,'nwave');
curWave = oiGet(oi,'wave');
meanIll = oiGet(oi,'meanilluminance');

if notDefined('newWave')
    prompt={'Start (nm)','Stop (nm)','Spacing (nm)'};
    def={num2str(curWave(1)),num2str(curWave(end)),num2str(oiGet(oi,'binwidth'))};
    dlgTitle='Wavelength resampling';
    lineNo=1;
    val =inputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(val), return; end
    
    low = str2num(val{1}); high = str2num(val{2}); skip = str2num(val{3});
    if high > low,       waveSpectrum.wave = low:skip:high;
    elseif high == low,  waveSpectrum.wave = low;     % User made monochrome, so onlyl 1 sample
    else 
        ieInWindowMessage('Bad wavelength ordering:  high < low. Data unchanged.',handles,5);
        return;
    end
else
    waveSpectrum.wave = newWave;
end

% if length(waveSpectrum.wave) > 1, 
%     waveSpectrum.binwidth = waveSpectrum.wave(2) - waveSpectrum.wave(1);
% else                    
%     waveSpectrum.binwidth = 1;
% end

% Current oi photons
photons = oiGet(oi,'photons');

% We clear the data to save memory space.  
oi = oiClearData(oi);

% We do this trick to be able to do a 1D interpolation. It is fast
% ... 2d is slow.  The RGB2XW format puts the photons in columns by
% wavelength.  The interp1 interpolates across wavelength
photons = RGB2XWFormat(photons)';
newPhotons = interp1(curWave,photons,waveSpectrum.wave)';
newPhotons = XW2RGBFormat(newPhotons,row,col);

oi = oiSet(oi,'spectrum',waveSpectrum); 
oi = oiSet(oi,'compressedphotons',newPhotons);

% Preserve the original mean luminance (stored in meanL) despite the resampling.
oi = oiSet(oi,'illuminance',oiCalculateIlluminance(oi));
oi = oiAdjustIlluminance(oi,meanL);

end