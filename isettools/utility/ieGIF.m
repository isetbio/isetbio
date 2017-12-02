function ieGIF(data, varargin)
% Save movie data (x, y, t) as a GIF
%
% Syntax:
%   [data, vObj] = ieGIF(data, [step], [vname], [delay], [hf])
% 
% Description:
%    Save the movie contained in data as a movie into a GIF
%
% Inputs:
%    data  - (row, col, color, time) or (row, col, time) (Required)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    step  - How many times frames to step over. Default = 1
%    vname - (video file name) Default = '/local/test.gif'
%    Delay - (pause between frames) Default = 0.03
%    hf    - Figure for showing data Default = vcNewGraphWin()
%

% History:
%    xx/xx/16 bw, jrg  ISETBIO Team (BW/JRG) 2016
%    11/22/17  jnm  Formatting. Added support for Windows in call - need to
%    check if this is necessary?
%

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

if isunix || ismac
    p.addParameter('vname', [isetbioRootPath '/local/test.gif'], @ischar);
else
    p.addParameter('vname', [isetbioRootPath '\local\test.gif'], @ischar);
end

p.parse(data, varargin{:});
data  = p.Results.data;
step  = p.Results.step;
vname      = p.Results.vname;
delay  = p.Results.delay;

%% Make GIF
figure(1)
for fr1 = 1:step:size(data, 3)
    
    imnew(:, :, 1, :) = data(:, :, fr1);
    imagesc(imnew)
    colormap gray;
    axis image;
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
      
    [imind, cm] = rgb2ind(im, 256);
    if fr1 == 1
        imwrite(imind, cm, vname, 'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, vname, 'gif', 'WriteMode', 'append', ...
            'DelayTime', delay);
    end
end