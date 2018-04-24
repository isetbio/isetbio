function [freq, fData] = ieSpace2Amp(pos, data, scaleData)
% Transform spatial data to amplitudes in cycles per spatial unit  
%
% Syntax:
%   [freq, fData] = ieSpace2Amp(pos, data, scaleData)
%
% Description:
%    The returned amplitudes are abs(fft(ifftshift(data))) of the data
%    vector. The spatial frequency value units are in cycles per unit of
%    the input data. For example, if the input data are in mm, then the
%    output data are in cycles per mm. If the input data are in meters,
%    then the output is cycles per meter.
%
%    The use of the ifftshift is because of the convention that these are
%    data as a function of space, with the natural position of the origin
%    being in the center of the input vector.  Because we are returning
%    only the amplitude, the ifftshift has no effect. But it would matter
%    if we were looking at the phase, too, so it seems wise to include it
%    as a matter of convention.
%
% Inputs:
%    pos       - Positions in spatial units (e.g., microns)
%    data      - vector of data values for each spatial position
%    scaleData - (Optional) Whether or not to scale the result values. If
%                this is true, the amplitude is normalized to a peak of 1.
%                Default is false.
%
% Outputs:
%    freq      - Frequency
%    fData     - Frequency data
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    scenePlot, oiPlot, sensorPlotLine
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/21/17  jnm  Formatting
%    12/30/17  dhb  Add in ifftshift.
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    data = sceneGet(scene, 'luminance');
    lum = data(5, :);
    pos = sceneSpatialSupport(scene, 'mm');
    [freq, fftlum] = ieSpace2Amp(pos.x, lum, 1);
    vcNewGraphWin;
    plot(freq,fftlum);
    grid on
    xlabel('Frequency (cy/mm)');
    ylabel('Amp');
%}

if notDefined('pos'), errordlg('You must define positions'); end
if notDefined('data'), errordlg('You must define a vector of data'); end
if notDefined('scaleData'), scaleData = 0; end

nSamp = length(data);
fData = abs(fft(ifftshift(data)));

% Uncomment this code to persuade yourself that the ifftshift has no
% effect.
% tData = abs(fft(data));
% if (max(abs(tData(:) - fData(:))) > 1e-10)
%     error('Adding in ifftshift has unexpected effect');
% end

% Scale the data to a peak of 1 before analyzing. Was default. Uh oh.
% Keep an eye on what has changed.
if scaleData, fData = fData / max(fData(:)); end

% This is the units per image 
unitPerImage = (max(pos) - min(pos));

% The frequency values are scaled into units. Without scaling the frequency
% representation is cycles/image (or cycle/data set). With this
% correction we have
%    cycle/image/(unit/image) = cycle/unit
% In addition, the true frequency numbers range from 0 (mean) on up, and
% there are only half as many (Nyquist) as there are samples.
freq = ((1:nSamp) - 1) / unitPerImage;
nFreq = round((nSamp - 1) / 2);

% Here are the frequency data from the mean up to the Nyquist frequency.
fData = fData(1:nFreq);
freq = freq(1:nFreq);

end