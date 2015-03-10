function [hdl, filterRGB] = sensorShowCFA(sensor,fullArray)
%Create an image illustrating the sensor CFA spatial pattern
%
%    [hdl, filterRGB] = sensorShowCFA(sensor,[fullArray = 0])
%
% Typically an image of just the pattern is shown.  The fullArray flag
% makes an image showing the full pattern.
%
% Example:
%   sensor = sensorCreate;
%   sensor = sensorSet(sensor,'pattern',[1 2 3; 2 3 1; 3 1 2]);
%   hdl = sensorShowCFA(sensor);
%
% Copyright ImagEval Consultants, LLC, 2010

if notDefined('sensor'),    sensor = vcGetObject('sensor'); end
if notDefined('fullArray'), fullArray = 0; end

% Get sensor information
pattern       = sensorGet(sensor,'pattern');
filterSpectra = sensorGet(sensor,'filterSpectra');
wave          = sensorGet(sensor,'wave');

if fullArray
    % Basically a call to sensorImageColorArray

    [cfaImage,mp] = sensorImageColorArray(sensorDetermineCFA(sensor));

    % This isn't quite as pretty for the human case.  
    %     vcNewGraphWin;
    %     image(cfaImage); colormap(mp)
    s = 3;
    filterRGB = imageIncreaseImageRGBSize(cfaImage,s);
    
    % Draw the CFA in a new figure window
    hdl = vcNewGraphWin; 
    set(hdl,'Name', 'CFA','menubar','None');
    image(filterRGB), colormap(mp); axis image off; truesize(hdl)
    
else
    % The human case usually drops into here.  The pattern and the full
    % array are really the same in that case, because for the human cone
    % mosaic there is no repeating block.

    % Create an RGB image illustrating the CFA colors
    % Possibly, we should be using the filter names rather than this
    [nRows,nCols] = size(pattern);
    filterRGB = zeros(nRows,nCols,3);
    irExtrapolation = 0.2;
    bMatrix = colorBlockMatrix(wave,irExtrapolation);
    for ii = 1:nRows
        for jj = 1:nCols
            colorFilter = filterSpectra(:, pattern(ii,jj));
            rgb = bMatrix'*colorFilter;
            if max(rgb) < 0.15
                % Probably IR, so give it a dark reddish gray
                filterRGB(ii,jj,:) = [0.4 0.3 0.2];
            else
                filterRGB(ii,jj,:) = (rgb'/max(rgb(:)));
            end
        end
    end

    % If it is small (not a cone mosaic) make it bigger.  Say 128 x 128
    % Otherwise, just make it bigger so we can see the individual points.
    % There is a chance that making it bigger runs us out of memory.  Sigh.
    if max(size(filterRGB,1)) < 64, s = 128/round(size(filterRGB,1));
    else                            s = 3;
    end
    filterRGB = imageIncreaseImageRGBSize(filterRGB,s);

    % Draw the CFA in a new figure window
    hdl = vcNewGraphWin; set(hdl,'Name', 'CFA','menubar','None');
    imagesc(filterRGB), axis image off; truesize(hdl)
end

return
