function [data, vObj] = ieMovie(data,varargin)
% Show a movie of an (x,y,t) or (x,y,c,t) matrix
%
%   [data, vObj] = ieMovie(data,varargin)
%
%  data:   (row,col,color,time) or (row,col,time) (Required)
%
% Name-value pairs
%                                             Default
%  step:   How many times frames to step over. (1);
%  show:   Display the movie                   (true)
%  save:   Write video to a file
%  vname:  Video file name                     (vName)
%  FrameRate: Video frame rate                 (20) - Only for video object
%  hf:        Figure for showing data          (vcNewGraphWin())
%  gamma:     Gamma exponent (d.^gamma)
%  ax:        Display axis
%
% Example:
%  ieMovie(rand(50,50,20),'show',false);
%  ieMovie(rand(50,50,20),'show',true,'gamma',1/1.5);
%  ieMovie(rand(50,50,20),'show',true,'save',true);
%  ieMovie(rand(50,50,20),'show',true,'save',true,'FrameRate',5);
%
%  [mov,vObj] = ieMovie(randn(50,50,3,20));
%
% ISETBIO Team (BW) 2016


%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('data',@isnumeric);

% Name-Value
p.addParameter('save',false,@islogical);
p.addParameter('vname','videoName',@ischar);
p.addParameter('FrameRate',20,@isnumeric);
p.addParameter('step',1,@isnumeric);
p.addParameter('show',true,@islogical);
p.addParameter('gamma',1,@isnumeric);
p.addParameter('ax',gca,@isgraphics);

p.parse(data,varargin{:});
data  = p.Results.data;
step  = p.Results.step;
save  = p.Results.save;
vname      = p.Results.vname;
show       = p.Results.show;
FrameRate  = p.Results.FrameRate;
gamma      = p.Results.gamma;
ax         = p.Results.ax;

%% Create the movie and video object

vObj = [];    % Video object
axes(ax);     % Select the axes

% Could be monochrome or rgb
tDim = ndims(data);     % Time step is always the last dimension
nFrames = size(data, tDim);

% Scale and gamma correct mov
if gamma ~= 1, data = ieScale(data,0,1) .^ gamma;
else,          data = ieScale(data,0,1);
end

% mind = min(data(:)); maxd = max(data(:));
mind = 0; maxd = 1;

% Create the video object if we plan to save
if save
    vObj = VideoWriter(vname);
    vObj.FrameRate = FrameRate;
    open(vObj);
end

% Step through each frame, saving in the video object (or not) and shownig
% on the screen (or not)
if isequal(tDim,4)
    % RGB data
    for ii=1:step:size(data,tDim)
        imagesc(data(:,:,:,ii)); axis image; set(gca,'xticklabel','','yticklabel','');
        caxis([mind maxd]); drawnow;
        if ii == 1 && ~show, set(gcf,'Visible','off'); end
        if save,  F = getframe; writeVideo(vObj,F); end
        pause(1/FrameRate);
    end
elseif isequal(tDim,3)
    % Monochrome data
    colormap(gray(256)); 
    for ii = 1:nFrames
        imagesc(data(:,:,ii)); axis image; set(gca,'xticklabel','','yticklabel','');
        caxis([mind maxd]); drawnow;
        if ii == 1 && ~show, set(gcf,'Visible','off'); end
        if save,  F = getframe; writeVideo(vObj,F); end
        pause(1/FrameRate);
    end
end

% Write the video object if save is true
if save
    writeVideo(vObj,F);
    close(vObj);
end

end
