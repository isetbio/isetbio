function gifName = ieGIF(data, varargin)
% Save movie data (x, y, t) as a GIF
%
% Syntax:
%   fileName = ieGIF(data, [step], [vname], [delay], [hf])
% 
% Description:
%    Save the movie contained in data into a GIF
%
% Inputs:
%    data  - (row, col, color, time) or (row, col, time) (Required)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    step    - How many times frames to step over. (Default = 1)
%    gifName - gif file name  (Default = 'isetbioRootPath/local/test.gif')
%    Delay   - Delay between frames (Default = 0.03 sec)
%    hf      - Figure for showing data Default = vcNewGraphWin()
%

% History:
%    xx/xx/16  bw, jrg  ISETBIO Team (BW/JRG) 2016
%    11/22/17  jnm      Formatting. Added support for Windows in call -
%                       need to check if this is necessary?

% Example:
%{
    ieGIF(rand(50, 50, 50));
%}

%% Parse inputs
p = inputParser;
p.addRequired('data', @isnumeric);
p.addParameter('delay', .03, @isnumeric);
p.addParameter('step', 1, @isnumeric);
p.KeepUnmatched = true;

p.addParameter('gifName', fullfile(isetbioRootPath,'local','test.gif'), @ischar);

p.parse(data, varargin{:});
data    = p.Results.data;
step    = p.Results.step;
gifName = p.Results.gifName;
delay   = p.Results.delay;

%% Make GIF
vcNewGraphWin;
cm = gray(max(data(:)));
colormap(cm); axis image; axis off;

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

end