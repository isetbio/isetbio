function gifName = ieGIF(data, varargin)
% Save movie data (x, y, t) as a gray scale GIF
%
% Syntax:
%   fileName = ieGIF(data, ...)
% 
% Description:
%    Save the movie in data as a GIF. Only works for gray scale.
%
%    Examples are included within the code.
%
% Inputs (required):
%    data    - (row, col, time) (Required)
%
% Outputs:
%    gifName - The name of the GIF
%
% Optional key/value pairs:
%    delay   - timing between frames in sec (default 0.05 sec)
%    gifName - file name (extension must be .gif, 
%              default is fullfile(isetbioRootPath, 'local', 'test.gif'))
%

% History:
%    xx/xx/16  bw, jrg  ISETBIO Team (BW/JRG) 2016
%    11/22/17  jnm      Formatting. Added support for Windows in call -
%                       need to check if this is necessary?
%    12/31/17   BW      Rewrote
%    01/16/18  jnm      Formatting update to match Wiki.

% Example:
%{
    ieGIF(ceil(255*rand(50, 50, 50)), 'delay', 0.3);
%}

%{
https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
So it appears that `imwrite` now supports writing animated grayscale GIFs
"When writing multiframe GIF images, X should be an 4-dimensional
M-by-N-by-1-by-P array, where P is the number of frames to write."

But if I try to pass it an M-by-N-by-3-by-P it seems to treat each RGB
color channel as a separate grayscale frame. Is there now way to write an
animated color GIF without a for loop over the frames?
%}

%% Parse inputs
p = inputParser;
p.addRequired('data', @isnumeric);
p.addParameter('delay', .05, @isnumeric);
p.KeepUnmatched = true;

p.addParameter('gifName', fullfile(isetbioRootPath, 'local', 'test.gif'), ...
    @ischar);

p.parse(data, varargin{:});
data    = p.Results.data;
gifName = p.Results.gifName;
delay   = p.Results.delay;

%% Make GIF modern way
[r, c, w]=size(data);
data = reshape(data, r, c, 1, w);
imwrite(data, gifName, 'gif', 'Loopcount', inf, 'DelayTime', delay);

end

%% Make GIF - from older examples using colormap and rgb2ind
%{
    vcNewGraphWin;
    cm = gray(max(data(:)));
    colormap(cm);
    axis image;
    axis off;

    for fr1 = 1:step:size(data, 3)
        imnew(:, :, 1, :) = data(:, :, fr1);
        image(imnew);
        % imagesc(imnew)
        % colormap gray;
        axis image; drawnow;
        frame = getframe(1);
        im = frame2im(frame);

        [imind, cm] = rgb2ind(im, 256);
        if fr1 == 1
            imwrite(imind, cm, gifName, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, gifName, 'gif', 'WriteMode', 'append', ...
            'DelayTime', delay);
        end
    end
%}
