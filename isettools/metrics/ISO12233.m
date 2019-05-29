function [results, fitme, esf, h] = ...
    ISO12233(barImage, deltaX, weight, plotOptions)
% ISO 12233 (slanted bar) spatial frequency response (SFR) analysis.
%
% Syntax:
%   [results, fitme, esf, h] = ...
%       ISO12233(barImage, deltaX, weight, plotOptions);
%
% Description:
%    Slanted-edge and color mis-registration analysis.
%
%    Run like this: ISO12233;
%    (You are prompted for a bar file and other parameters)
%
%    This function contains example of usage. To access, enter 'edit
%
% Inputs:
%    barImage    - VARIES. Either a matrix representing the RGB image of
%                  the slanted bar, or a structure which contains the
%                  matrix image.
%    deltaX      - (Optional) Numeric. The sensor sample spacing in
%                  millimeters (expected). It is possible to send in a
%                  display spacing in dots per inch (dpi), in which case
%                  the number is > 1 (it never is for sensor sample
%                  spacing). In that case, the value returned is cpd on the
%                  display at a 1m viewing distance. Default 0.002.
%    weight      - (Optional) Matrix. A 1x3 matrix of the luminance
%                  weights; in RGB order. Default [0.3, 0.6, 0.1].
%    plotOptions - (Optional) String. A string from the set of {'all',
%                  'luminance', or 'none'}. Default 'all'.
%
% Outputs:
%    results     - Structure. The created results structure.
%    fitme       - Matrix. A matrix representing the full linear fit to
%                  something... (to what?)
%    esf         - Matrix. Appears to be an interpolated edge for each of
%                  the color channels.
%    h           - Figure. Appears to be the resulting figure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * General notes:
%      - One cycle on the sensor has a frequency of:
%        1 / sensorGet(sensor, 'width', 'mm')
%      - One cycle on the sensor is 1 / sensorGet(sensor, 'fov') cycles in
%        the original image.
%    * PROGRAMMING TODO: Rather than decode pixelWidth and dpi based on the
%      value, we should probably set a flag and be explicit.
%
% References:
%    This code originated with Peter Burns, peter.burns@kodak.com
%    12 August 2003, Copyright (c) International Imaging Industry
%    Association, and Substantially re-written by ImagEval Consulting, LLC
%

% Hisory:
%    08/12/03  PB   Original (Peter Burns - International Imaging Industry)
%    XX/XX/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/17/18  dhb  Deleted broken example that relied on sensorGet.
%                   Skip auto-execute of other examples, as they need
%                   user input and/or print warnings.

% Examples:
%{
    % ETTBSkip. This is interactive and should not be autorun. But some
    % comments about what user should do would be great.

    % Interactive usage 
    ISO12233;
%}
%{
    % ETTBSkip. This also requiers user input and should not be autorun.
    % But some comments about what user should do would be great.

    deltaX = 0.006;  % Six micron pixel. deltaX Units appear to be mm.
    [results, fitme, esf] = ISO12233([], deltaX, []);

    % The whole thing and cpd assuming a 1m viewing distance
    rectMTF = [xmin ymin width height];
    c = rectMTF(3) + 1;
    r = rectMTF(4) + 1;
    roiMTFLocs = ieRoi2Locs(rectMTF);
    barImage = vcGetROIData(vciBlurred, roiMTFLocs, 'results');
    barImage = reshape(barImage, r, c, 3);
    wgts = [ 0.3 0.6 0.1];
    [results, fitme, esf] = ISO12233(barImage, deltaX, wgts, 'luminance');
%}

%%
if notDefined('deltaX')
    deltaX = .002;
    warning('Assuming 2 micron pixel');
end
if notDefined('weight')
    % RGB: Luminance weights
    weight = [0.3, 0.6, 0.1];
end
if notDefined('plotOptions')
    % Options are: all, luminance, or none
    plotOptions = 'all';
end
if notDefined('barImage')
    % If there is no image, then you can read a file with the bar image.
    % You are also asked to specify a look-up table file that converts the
    % data in the edgeFile into linear units appropriate for the MTF
    % calculation. If no lutFile is selected, then we assume this
    % transformation is not necessary for your data.
    edgeFile = vcSelectDataFile('stayput', 'r', [], ...
        'Select slanted bar image');
    [barImage, smax] = readBarImage(edgeFile);
    lutFile = vcSelectDataFile('stayput', 'r', [], 'Select LUT file');

    if ~isempty(lutFile)
        [oepath, oename, oeext] = fileparts(lutFile);
        % Convert through LUT and make sure data are in double format
        barImage = getoecf(barImage, oepath, [oename, oeext]);
    end
    barImage = double(barImage);
elseif isstruct(barImage) && isequal(barImage.type, 'vcimage')
    % The barImage is really the image processor (ip)
    ip = barImage;
    rect = ISOFindSlantedBar(ip); 
    roiLocs = ieRoi2Locs(rect);
    barImage = vcGetROIData(ip, roiLocs, 'results');
    col = rect(3) + 1;
    row = rect(4) + 1;
    barImage = reshape(barImage, row, col, 3);
    smax = max(barImage(:));
else
    smax = max(barImage(:));
    % edgeFile = 'Input data';
end

% Extract region of interest
[nlow, nhigh, cstatus] = clipping(barImage, 0, smax, 0.005);
if cstatus ~= 1
    fprintf('Fraction low data: %.3f\n', nlow);
    fprintf('Fraction high data: %.3f\n', nhigh);
end

% Default sampling and color weights
if deltaX == 1
    % Unknown physical units. Not preferred.
    funit = 'cy/pixel';
elseif deltaX > 1
    % Spacing is with respect to the display RGB (dpi). We assume a 1m
    % viewing distance.
    deltaX = 25.4/deltaX;
    funit = 'cy/deg at 1m distance';
else
    % Spacing is with respect to the sensor surface. Value should be in mm.
    funit = 'cy/mm on sensor';
end

% Default:  weight = [0.3, 0.6, 0.1];
[nRow, nCol, nWave] = size(barImage);
if (nRow < 5) || (nCol < 5), warning('Image region too small'); return; end

%% Start computations
% Form luminance record using the weight vector for red, green and blue
% Add this as a fourth image to barImage
if nWave == 3
    % lum = zeros(nRow, nCol);
    lum = weight(1) * barImage(:, :, 1) + ...
        weight(2) * barImage(:, :, 2) + weight(3) * barImage(:, :, 3);
    barImage(:, :, 4) = lum;
    nWave = 4;
end

% rotate horizontal edge to vertical
[barImage, nRow, nCol, rflag] = rotatev(barImage);
loc = zeros(nWave, nRow);

% Need 'positive' edge for good centroid calculation
fil1 = [0.5 -0.5];
fil2 = [0.5 0 -0.5];
tleft = sum(sum(barImage(:, 1:5, 1), 2));
tright = sum(sum(barImage(:, nCol - 5:nCol, 1), 2));
if tleft > tright
    fil1 = [-0.5 0.5];
    fil2 = [-0.5 0 0.5];
end

% Test for low contrast edge;
test = abs((tleft - tright) / (tleft + tright));
if test < 0.2
    disp(' ** WARNING: Edge contrast is less that 20%, this can');
    disp('             lead to high error in the SFR measurement.');
end

fitme = zeros(nWave, 2);
slout = zeros(nWave, 1);

% smoothing window for first part of edge location estimation -
% to used on each line of ROI
win1 = ahamming(nCol, (nCol+1)/2);      % Symmetric window
for color = 1:nWave                       % Loop for each color
    %     if nWave == 1, pname = ' ';
    %     else pname = [' Red ' 'Green'  'Blue ' ' Lum '];
    %     end
    c = deriv1(barImage(:, :, color), nRow, nCol, fil1);
    % vcNewGraphWin;
    % imagesc(c);
    % colormap(gray)
    % compute centroid for derivative array for each line in ROI.
    % NOTE WINDOW array 'win'
    for n = 1:nRow
        % -0.5 shift for FIR phase
        loc(color, n) = centroid(c(n, 1:nCol)' .* win1) - 0.5;
    end
    % clear c

    fitme(color, :) = findedge(loc(color, :), nRow);
    place = zeros(nRow, 1);
    for n = 1:nRow
        place(n) = fitme(color, 2) + fitme(color, 1) * n;
        win2 = ahamming(nCol, place(n));
        loc(color, n) = centroid(c(n, 1:nCol)' .* win2) -0.5;
    end

    fitme(color, :) = findedge(loc(color, :), nRow);
    % fitme(color, :); % used previously to list fit equations
end

summary{1} = ' ';  % initialize
nWaveOut = nWave;  % output edge location listing
if nWave == 4, nWaveOut = nWave - 1; end

midloc = zeros(nWaveOut, 1);
summary{1} = 'Edge location, slope'; % initialize

for i = 1:nWaveOut
    % slope is as normally defined in image coords.
    slout(i) = - 1 ./ fitme(i, 1);
    % positive flag it ROI was rotated
    if rflag == 1, slout(i) = -fitme(i, 1); end

    % evaluate equation(s) at the middle line as edge location
    midloc(i) = fitme(i, 2) + fitme(i, 1) * ((nRow - 1) / 2);

    summary{i + 1} = [midloc(i), slout(i)];
end

% Could insert a display flag
% disp('Edge location(s) and slopes = ');
% disp([midloc(1:nWaveOut), slout(1:nWaveOut)]);
if nWave>2
    summary{1} = strcat('Edge location, slope, misregistration ', ...
        '(second record, G, is reference)');
    misreg = zeros(nWaveOut, 1);
    for i = 1:nWaveOut
        misreg(i) = midloc(i) - midloc(2);
        summary{i + 1} = [midloc(i), slout(i), misreg(i)];
    end

    % Turned off display
    % disp('Misregistration, with green as reference (R, G, B, Lum) = ');
    % for i = 1:nWaveOut
    %     fprintf('%10.4f\n', misreg(i));
    % end
end

% Full linear fit is available as variable fitme. Note that the fit is for
% the projection onto the X-axis, 
%       x = fitme(color, 1) y + fitme(color, 2)
% so the slope is the inverse of the one that you may expect
nbin = 4;
nn = floor(nCol * nbin);
mtf = zeros(nn, nWave);
nn2 = nn / 2 + 1;

% We compute the frequencies in terms of the deltaX spacing. Ordinarily, 
% deltaX is in terms of the pixel pitch on the sensor in millimeters. Some
% convenient facts:
%    1 cycle on the sensor is the frequency:
%        1 / sensorGet(sensor, 'width', 'mm')
%    1 cycle on the sensor is 1 / sensorGet(sensor, 'fov') cycles in the
%        original image.
%    1 cycle in the barImage is 1/barImageWidthDeg
%
% If the deltaX input is > 1, however, we interpret the parameter as
% display dpi. In this case the cycles are converted from cy/mm to cy/deg
% assuming a viewer located 1 m from the display.
%
% The first step in the code converts the frequencies into lines per
% millimeter on the sensor or display surface.
freq = zeros(nn, 1);
for n = 1:nn, freq(n) = nbin * (n - 1) / (deltaX * nn); end
% limits plotted sfr to 0-1 cy/pxel freqlim = 2 for all data
freqlim = 1;
nn2out = round(nn2 * freqlim / 2);  % What is this?
nfreq = n / (2 * deltaX * nn);      % half-sampling (Nyquist) frequency

% If the units are for the display surface, we further convert to cycles
% per degree on the display by converting lpm and assuming a viewing
% distance of 1M
if strncmpi(funit, 'cy/deg', 6)
    % We currently have cyc/mm, so we need mm/deg for the conversion.
    %
    % The right angle from the viewer to screen produces the equation:
    % Opposite = X; Adjacent = 1;
    % tan(X/adjacent) = ieDeg2rad(0.5)
    % X = 2 * atan(ieDeg2rad(0.5)) * adjacent  is the answer in meters.
    % Multiply by 1000 to get mm per deg on the display surface
    mmPerDeg = 2 * atan(ieDeg2rad(0.5)) * 1e3;

    % cyc/mm * mm/deg give us cyc/deg
    freq = freq * mmPerDeg;
    nfreq = nfreq * mmPerDeg;
end

% centered Hamming window
win = ahamming(nbin * nCol, (nbin * nCol + 1) / 2);

% Loop for each color record. This variable is returned. It seems to be the
% interpolated edge for each of the color channels.
esf = zeros(nn, nWave);

for color = 1:nWave
    % project and bin data in 4x sampled array
    point = project(barImage(:, :, color), loc(color, 1), ...
        fitme(color, 1), nbin);
    % vcNewGraphWin; plot(point); colormap(gray)
    esf(:, color) = point;  % Not sure what esf stands for.
                            % Estimated spatial frequency?

    % compute first derivative via FIR (1x3) filter fil
    c = deriv1(point', 1, nn, fil2);  % vcNewGraphWin; plot(c)
    c = c';
    mid = centroid(c);
    temp = cent(c, round(mid));       % shift array so it is centered
    c = temp;
    clear temp;

    % apply window (symmetric Hamming)
    c = win .* c;    % vcNewGraphWin; plot(c)

    % Transform, scale %% The FFT of the point spread, 
    temp = abs(fft(c, nn));    % vcNewGraphWin; plot(temp)
    mtf(1:nn2, color) = temp(1:nn2) / temp(1);
end

dat = zeros(nn2out, nWave + 1);
for i = 1:nn2, dat(i, :) = [freq(i), mtf(i, :)]; end

%% Plot SFRs on same axes
if nWave > 1
    sym{1} = '-r';
    sym{2} = '-g';
    sym{3} = '-b';
    sym{4} = '-k';
else
    sym{1} = 'k';
end
% ttext = sprintf('ISO 12233: %s', edgeFile);

% screen = get(0, 'ScreenSize');
% % defpos = get(0, 'DefaultFigurePosition');
% set(0, 'DefaultFigurePosition', [15 25 0.6*screen(3) 0.4*screen(4)]);

results.freq = freq(1:nn2out);
results.mtf = mtf(1:nn2out, :);
results.nyquistf = nfreq;
lumMTF = results.mtf(:, end);

belowNyquist = (results.freq < nfreq);

% Sometimes, if the image is very noisy, lumMTF has a number of NaNs. We
% won't find mtf50 in such cases.
if ~isnan(lumMTF)
    % Old calculation
    %   results.mtf50 = interp1(lumMTF, results.freq, 0.5);
    % New calculation
    %   Sample freq finely
    %   Find the below-nyquist freq closest to an MTF value of 0.5
    iFreq = 0:0.2:results.nyquistf;
    iLumMTF = interp1(results.freq(belowNyquist), ...
        lumMTF(belowNyquist), iFreq);
    [v, idx] = min(abs(iLumMTF - 0.5));
    results.mtf50 = iFreq(idx);
else
    fprintf('NaN lumMTF values. No plot is generated.\n');
    return
end

% The area under the curve to the right of the nyquist as a percentage of
% the total area in the green channel (when RGB), or in the luminance
% (first channel) when a monochrome image.
if nWave == 4
    results.aliasingPercentage = 100 * ...
        sum(results.mtf(~belowNyquist, 2)) / sum(results.mtf(:, 2));
elseif nWave == 1
    results.aliasingPercentage = 100 * ...
        sum(results.mtf(~belowNyquist, 1)) / sum(results.mtf(:, 1));
end

switch plotOptions
    case 'all'
        % Set data into the figure
        h = vcNewGraphWin;
        set(h, 'userdata', results);
        % Draw the luminance term
        p = plot(freq(1:nn2out), mtf(1:nn2out, 1), sym{1});
        set(p, 'linewidth', 2);

        title('ISO 12233');
        xlabel(['Spatial frequency (', funit, ')']);
        ylabel('Contrast reduction (SFR)');

        hold on;
        if nWave > 1
            for n = 2:nWave
                p = plot(freq(1:nn2out), mtf(1:nn2out, n), sym{n});
                set(p, 'linewidth', 2);
            end
        end

        % Half-sampling line on graph
        line([nfreq , nfreq], [0.05, 0]), 

        % TODO: Put little box or legend, or some kind lines/points to
        % indicate these on the graph.
        txt1 = sprintf('Nyquist = %0.2f\n', nfreq);
        txt2 = sprintf('Mtf50 = %0.2f\n', results.mtf50);
        txt3 = sprintf('Percent alias = %0.2f\n', ...
            results.aliasingPercentage);
        txt = addText(txt1, txt2);
        txt = addText(txt, txt3);

        % delta is (x, y)
        plotTextString(txt, 'ur', [0.4 0.2], 18);
        hold off;
        grid on;
    case 'luminance'
        % Set data into the figure
        h = vcNewGraphWin;
        set(h, 'userdata', results);

        p = plot(freq(1:nn2out), mtf(1:nn2out, 1), sym{1});
        set(p, 'linewidth', 2);

        % [p, fname, e] = fileparts(ttext);
        title('ISO 12233');
        xlabel(['Spatial frequency (', funit, ')']);
        ylabel('Contrast reduction (SFR)');

        % Half-sampling line on graph
        line([nfreq , nfreq], [0.05, 0]), 

        % TODO: Put little box or legend, or some kind lines/points to
        % indicate these on the graph.
        % See above for a fix to this.
        txt = sprintf('Nyquist = %0.2f\n', nfreq);
        newText = sprintf('Mtf50 = %0.2f\n', results.mtf50);
        txt = addText(txt, newText);
        newText = sprintf('Percent alias = %0.2f', ...
            results.aliasingPercentage);
        txt = addText(txt, [newText, ' %']);
        plotTextString(txt, 'ur');
        hold off;
    case 'none'
        % Do nothing, we don't want a plot
        h = [];
    otherwise
        error('Unknown plotOptions: %s\n', plotOptions);
end
% return;
end

%---------------------------------------------------
% ISO 12233 subroutines
%---------------------------------------------------
function [barImage, smax] = readBarImage(edgeFile)
% Read the image file.
%
% Syntax:
%   [barImage, smax] = readBarImage(edgeFile)
%
% Description:
%    Read and return the bar image contained in edgeFile.
%
% Inputs:
%    edgeFile - String. A string containing the filepath for the bar image.
%
% Outputs:
%    barImage - Matrix. The barImage matrix.
%    smax     - Numeric. The maximum size of the image based on image type.
%
% Optional key/value pairs:
%    None.
%

tempBarImage = imread(edgeFile);
% [nRow nCol nWave] = size(tempBarImage);

switch (lower(class(tempBarImage)))
    case 'uint8'
        smax = 255;
    case 'uint16'
        smax = 2 ^ 16 - 1;
    otherwise
        smax = 1e10;
end
barImage = getroi(tempBarImage);
% return;
end

%----------------------------------------------------
function [data] = ahamming(n, mid)
% Generate an asymmetrical Hamming Window array.
%
% Syntax:
%   data = ahamming(n, mid)
%
% Description:
%    [data] = ahamming(n, mid)  Generates asymmetrical Hamming window
%    array. If mid = (n + 1) / 2 then the usual symmetrical Hamming array
%    is returned.
%
% Inputs:
%    n    - Numeric. The length of the array.
%    mid  - Numeric. The midpoint (maximum) of the window function.
%
% Outputs:
%    data - Matrix. An nx1 window array.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

data = zeros(n, 1);

wid1 = mid - 1;
wid2 = n - mid;
wid = max(wid1, wid2);
for i = 1:n
    arg = i - mid;
    data(i) = 0.54 + 0.46 * cos(pi * arg / wid);
end
% return;
end

%----------------------------------------------------
function [b] = cent(a, center)
% A one dimensional array shift for centering data
%
% Syntax:
%   b = cent(a, center)
%
% Description:
%    Shift one-dimensional array, so that a(center) is located at
%    b(round((n + 1) / 2).
%    Written to shift a line-spread function array prior to
%    applying a smoothing window.
%
% Inputs:
%    a      - Array. The Input array.
%    center - Numeric. The signal center to be shifted to.
%
% Outputs:
%    b      - Array. The shifted array.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

n = length(a);
b = zeros(n, 1);
mid = round((n + 1) / 2);

del = round(center - mid);

if del > 0
    for i = 1:n - del, b(i) = a(i + del); end
elseif del < 1
    for i = -del + 1:n, b(i) = a(i + del); end
else
    b = a;
end

end

%----------------------------------------------------
function [loc] = centroid(x)
% Find the centroid of a provided vector.
%
% Syntax:
%   [loc] = centroid(x)
%
% Description:
%    Return the centroid location of a vector.
%
% Inputs:
%    x   - Array. A vector array.
%
% Outputs:
%    loc - Numeric. The centroid in units of array index.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

loc = 0;
for n = 1:length(x), loc = loc + n * x(n); end

if sum(x) == 0, warndlg('Values are all zero. Invalid centroid'); end
loc = loc/sum(x);
% return;
end

%----------------------------------------------------
function [nlow, nhigh, status] = clipping(barImage, low, high, thresh1)
% Check the data for clipping
%
% Syntax:
%   [nlow, nhigh, status] = clipping(a, low, high, thresh1)
%
% Description:
%    Function checks for clipping of data array.
%
% Inputs:
%    barImage - Vector. A vector representing the bar image.
%    low      - Numeric. The low clip value.
%    high     - Numeric. The high clip value.
%    thresh1  - Numeric. The threshold fraction [0-1] used for warning. If
%               set to 0, all clipping is reported.
%
% Outputs:
%    nlow     - Numeric. The number of low clippings.
%    nhigh    - Numeric. The number of high clippings.
%    status   - Boolean. A numeric boolean to represent whether or not
%               clipping errors have occurred.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

status = 1;

[nRow, nCol, nW] = size(barImage);
nhigh = zeros(nW, 1);
nlow = zeros(nW, 1);

for k = 1: nW
    for j = 1: nCol
        for i = 1: nRow
            if barImage(i, j, k) < low
                nlow(k) = nlow(k) + 1;
            end
            if barImage(i, j, k) > high
                nhigh(k) = nhigh(k) + 1;
            end
        end
    end
end

nhigh = nhigh ./ (nRow * nCol);
for k = 1: nW
    if nlow(k) > thresh1
        disp([' *** Warning: low clipping in record ', num2str(k)]);
        status = 0;
    end
    if nhigh(k) > thresh1
        disp([' *** Warning: high clipping in record ', num2str(k)]);
        status = 0;
    end
end
nlow = nlow ./ (nRow * nCol);
if status ~= 1, warndlg('Data clipping errors detected', 'ClipCheck'); end
% return;
end

%----------------------------------------------------
function  [b] = deriv1(a, nRow, nCol, fil)
% First derivative of the provided array
%
% Syntax:
%   b = deriv1(a, nRow, nCol, fil)
%
% Description:
%    This function computes first derivative via FIR (1xn) filter. It also
%    suppresses the edge effects and preserves the vector size. The filter
%    is applied in the nCol direction only.
%
% Inputs:
%    a    - Matrix. A data array matrix of nRow by nCol.
%    nRow - Numeric. The number of rows in a.
%    nCol - Numeric. The number of columns in a.
%    fil  - The array of filter coefficients (ex. [-0.5 0.5])
%
% Outputs:
%    b    - Matrix. A data array of nRow by nCol containing the filtered a.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

b = zeros(nRow, nCol);
nn = length(fil);
for i = 1:nRow
    temp = conv(fil, a(i, :));
    b(i, nn:nCol) = temp(nn:nCol);  %ignore edge effects, preserve size
    b(i, nn-1) = b(i, nn);
end
% return;
end

%----------------------------------------------------
function  [slope, int] = findedge(cent, nRow)
% Fit a linear equation to the data.
%
% Syntax:
%   [slope, int] = findedge(cent, nRow)
%
% Description:
%    Fit linear equation to data, written to process edge location array.
%
%    Using the least-square format to fit the line, the returns follow the
%    format of x = int + slope * cent(x). Please note that this is the
%    inverse of the typical cent(x) = int + slope * x calculation.
%
% Inputs:
%    cent  - Array. An array of the centroid values.
%    nRow  - Numeric. The length of the provided cent array.
%
% Outputs:
%    slope - Numeric. The slope of the least-square fit line.
%    int   - Numeric. The numeric offset from the axes for the least-square
%            fit line.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass.

index = 0:nRow - 1;
[slope, int] = polyfit(index, cent, 1);  % x = f(y)
% return
end

%----------------------------------------------------
function [array, status] = getoecf(array, oepath, oename)
% Read and apply the oecf
%
% Syntax:
%   [array, status] = getoecf(array, oepath, oename)
%
% Description:
%    Reads look-up table and applies it to a data array.
%
% Inputs:
%    array  - Matrix. A data array matrix of nRow by pnix by nWave
%    oepath - String. A filepath string for the oecf table.
%    oename - String. The filename containing the oecf table. Please note
%    that this is a tab-delimited text file for the table (256x1, 256, 3).
%
% Outputs:
%    array  - Matrix. The transformed input array.
%    status - Boolean. A boolean indicating whether or not there were
%             errors in applying the oecf. A status of 0 means no error. A
%             1 indicates a bad table file.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

status = 0;
stuff = size(array);
nRow = stuff(1);
nCol = stuff(2);
if isequal(size(stuff), [1 2]), nWave = 1; else, nWave = stuff(3); end

temp = fullfile(oepath, oename);
oedat = load(temp);
%oedat = oename;
dimo = size(oedat);
if dimo(2) ~= nWave, status = 1; return; end
if nWave == 1
    for i = 1:nRow
        for j = 1:nCol, array(i, j) = oedat(array(i, j) + 1, nWave); end
    end
else
    for i = 1:nRow
        for j = 1:nCol
            for k = 1:nWave
                array(i, j, k) = oedat(array(i, j, k) + 1, k);
            end
        end
    end
end
% return;
end

%----------------------------------------------------
function [select, coord] = getroi(array)
% Select and return a region of interest
%
% Syntax:
%   [select, coord] = getroi(array)
%
% Description:
%    Select and return a region of interest (ROI) from an image via a GUI
%    window and 'right-button-mouse' operation. If the mouse button is
%    clicked and released without movement, it will be treated as if the
%    entire displayed image was selected.
%
% Inputs:
%    array  - Matrix. An input image array of (nRow, nCol, [ncolor]), of
%             uint8 format. (ncolor is denoted as optional by brackets).
%
% Outputs:
%    select - Matrix. The output ROI as an array of the format (newlin,
%             newpix, [ncolor]), following the format specified by array.
%    coord  - Matrix. A list of coordinates of the ROI of the format
%             (upperleft(x, y), lowerright(x, y)).
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/23/19  JNM  Documentation pass

dim = size(array);
nRow = dim(1);
nCol = dim(2);
if isequal(size(dim), [1 2]), nWave = 1; else, nWave = dim(3); end
% 0.95 is to allow a tolerance which I need so very large narrow images
% stay visible
screen = 0.95 * (get(0, 'ScreenSize'));
% screen = get(0, 'ScreenSize');

% Set aspect ratio approx to that of array
rat = nCol / nRow;
% following lines make ROI selection of narrow images easier
if rat < 0.25, rat = 0.25; elseif rat > 4, rat = 4; end
if nRow >= nCol
    if nRow > 0.5 * screen(4)
        % This change helps with large images
        ht = min(nRow, 0.8 * screen(4));
    else
        ht = 0.5 * screen(4);
    end
    wid = round(ht * rat);
else
    if nCol > 0.5 * screen(3)
        % This change helps with large images
        wid = min(nCol, 0.8 * screen(3));
    else
        wid = 0.5 * screen(3);
    end
    ht = round(wid / rat);
end
pos = round([screen(3)/10 screen(4)/10 wid ht]);
figure(1), set(gcf, 'Position', pos);

disp(' ');
disp('Select ROI with right mouse button, no move = all');
disp(' ');

temp = class(array);
if ~strcmp(temp(1:5), 'uint8')
    imagesc(double(array) / double(max(max(max(array))))), 
    colormap('gray'), 
    title('Select ROI');
else
    if nWave == 1
        imagesc(array), 
        colormap('gray'), 
        title('Select ROI');
    else
        imagesc(array), 
        colormap('gray'), 
        title('Select ROI');
    end
end
%axis off

% junk = waitforbuttonpress;
waitforbuttonpress;
ul = get(gca, 'CurrentPoint');
% final_rect = rbbox;
rbbox;
lr = get(gca, 'CurrentPoint');
ul = round(ul(1, 1:2));
lr = round(lr(1, 1:2));

if ul(1, 1) > lr(1, 1)  % sort x coordinates
    mtemp = ul(1, 1);
    ul(1, 1) = lr(1, 1);
    lr(1, 1) = mtemp;
end
if ul(1, 2) > lr(1, 2)  % sort y coordinates
    mtemp = ul(1, 2);
    ul(1, 2) = lr(1, 2);
    lr(1, 2) = mtemp;
end

% if del x, y <10 pixels, select whole array
roi = [lr(2) - ul(2)  lr(1) - ul(1)];
if roi(1) < 10
    ul(2) = 1;
    lr(2) = nRow;
end
if roi(2) < 10
    ul(1) = 1;
    lr(1) = nCol;
end
select = double(array(ul(2):lr(2), ul(1):lr(1), :));
coord = [ul(:, :), lr(:, :)];
close;
% return;
end

%----------------------------------------------------
function [point, status] = project(barImage, loc, slope, fac)
% Projects data along the slanted edge to a common line
%
% Syntax:
%   [point, status] = project(barImage, loc, slope, fac)
%
% Description:
%    The data in array barImage are projected along the direction defined
%    by the function:
%        nCol = (1 / slope) * nRow
%    Data are accumulated in 'bins' that have a width (1 / fac) pixel.
%
%    The supersampled one-dimensional vector is returned.
%
% Inputs:
%    barImage - Matrix. The input data array of image
%    loc      - Array. An unused variable???
%    slope    - Numeric. The slope is calculated from the least-square fit
%               to the edge in a separate routine called findedge.
%    fac      - (Optional) Numeric. The factor. Default 4.
%
%    x = loc + slope*cent(x)
%         Note that this is the inverse of the usual cent(x) = int + slope*x
%
%  fac:    oversampling (binning) factor  (default = 4)
%  point:  output vector of
%  status = 1, OK
%  status = 0, zero counts encountered in binning operation
%              arning is printed, but execution continues
%
% See also: sfrmat11 and sfrmat2 functions.
%
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association
%
% Edited by ImagEval, 2006-2007
%
% Notes:  This so-called projection operator is really an alignment
% operation. The data along the different rows are slid so that the
% edge of the slanted bar on the different lines is the same. I
% think.

if (nargin < 4), fac = 4; end
status = 1;                      % Assume we are good to go
[nRow, nCol] = size(barImage);   % Size of the bar image
% figure(1); imagesc(barImage); axis image; colormap(gray(255))

% big = 0;
nn = nCol * fac;

% smoothing window. Why is this not used?  I commented out.
% win = ahamming(nn, fac*loc(1, 1));
% plot(win)
slope = 1 / slope;
offset = round(fac * (0  - (nRow - 1) / slope));

del = abs(offset);
if offset > 0, offset = 0; end
barray = zeros(2, nn + del + 100);

% Projection and binning
for n = 1:nCol
    for m = 1:nRow
        x = n - 1;
        y = m - 1;
        ling = ceil((x - y / slope) * fac) + 1 - offset;
        barray(1, ling) = barray(1, ling) + 1;
        barray(2, ling) = barray(2, ling) + barImage(m, n);
    end
end

% Initialize
point = zeros(nn, 1);
start = 1 + round(0.5 * del);

% Check for zero counts
nz = 0;
for i = start:start + nn - 1
    if barray(1, i) == 0
        nz = nz + 1;
        status = 0;
        if i == 1, barray(1, i) = barray(1, i + 1);
        else, barray(1, i) = (barray(1, i - 1) + barray(1, i + 1)) / 2;
        end
    end
end

if status == 0
    disp('                            WARNING');
    disp('      Zero count(s) found during projection binning. The edge ')
    disp('      angle may be large, you may need more lines of data.');
    disp('      Or the short edges of the rect may not cover the line');
    disp('      Execution will continue, but see Users Guide for info.');
    disp(nz);
end

% Combine into a single edge profile, point
for i = 0:nn - 1
    point(i + 1) = barray(2, i + start) / barray(1, i + start);
end

% This is the returned unified edge profile
% figure(1);
% plot(point);

% return;
end


%----------------------------------------------------
function [a, nRow, nCol, rflag] = rotatev(a)
% Rotate the provided array.
%
% Syntax:
%   [a, nRow, nCol, rflag] = rotatev(a)
%
% Description:
%    Rotate array so that long dimensions is vertical (line) drection.
%
% Inputs:
%    a     - Matrix. The matrix array to rotate.
%
% Outputs:
%    a     - Matrix. The modified input matrix.
%    nRow  - Numeric. The number of rows in the modified matrix.
%    nCol  - Numeric. The number of columns in the modified matrix.
%    rflag - Boolean. A boolean indicating whether a rotation happened. A 1
%            means that a rotation happened.
%
% Optional key/value pairs:
%    None.
%

% History:
%    08/05/02  PB   Peter Burns 5 Aug. 2002
%                   Copyright(c) International Imaging Industry Association
%    05/28/19  JNM  Documentation pass

dim = size(a);
nRow = dim(1);
nCol = dim(2);

if isequal(size(dim), [1 2]), nWave = 1; else, nWave = dim(3); end

rflag = 0;
if nCol > nRow
    rflag = 1;
    b = zeros(nCol, nRow, nWave);
    % temp = zeros(nCol, nRow);
    for i = 1:nWave
        temp = a(:, :, i)';
        b(:, :, i) = temp;
    end
    a = b;
    temp = nRow;
    nRow = nCol;
    nCol = temp;
end
% return
end
