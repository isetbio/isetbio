function ieGIF(data,varargin)
% Save movie data (x,y,t) as a GIF
%
%   [data, vObj] = ieGIF(data,varargin)
% 
%  data:   (row,col,color,time) or (row,col,time) (Required)
%  step:   How many times frames to step over. Default = 1;
%  vname:  (video file name)
%  Delay:  (pause between frames)
%  hf:        Figure for showing data (vcNewGraphWin() by default)
%
% Example:
%   ieGif(rand(50,50,50));
%
% ISETBIO Team (BW/JRG) 2016


%% Parse inputs
p = inputParser;
p.addRequired('data',@isnumeric);
p.addParameter('vname',[isetbioRootPath '/local/test.gif'],@ischar);
p.addParameter('delay',.03,@isnumeric);
p.addParameter('step',1,@isnumeric);
p.KeepUnmatched = true;

p.parse(data,varargin{:});
data  = p.Results.data;
step  = p.Results.step;
vname      = p.Results.vname;
delay  = p.Results.delay;

%% Make GIF
figure(1)
for fr1 = 1:step:size(data,3)
    
    imnew(:,:,1,:) = data(:,:,fr1);
    imagesc(imnew); colormap gray; axis image; drawnow;
    frame = getframe(1);
    
      im = frame2im(frame);
      
      [imind,cm] = rgb2ind(im,256);
    if fr1 == 1;
        imwrite(imind,cm,vname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,vname,'gif','WriteMode','append','DelayTime',delay);
    end
end