function [data, vObj] = ieMovie(data,varargin)
% Show a movie of an (x,y,t) or (x,y,c,t) matrix
%
%   [data, vObj] = ieMovie(data,varargin)
% 
%  data:   (row,col,color,time) or (row,col,time) (Required)
%  step:   How many times frames to step over. Default = 1;
%  show:   Display the movie
%  vname:  (video file name)
%  FrameRate: (video frame rate)
%  hf:        Figure for showing data (vcNewGraphWin() by default)
%
% Example:
%   ieMovie(rand(50,50,50));
%   
%   dFlag = true;
%   ieMovie(rand(50,50,50),'dFlag',dFlag);
%
%   dFlag.vname = 'tmp'; dFlag.FrameRate = 5; dFlag.hf = vcNewGraphWin;
%   [mov,vObj] = ieMovie(rand(50,50,50),'step',3,'dFlag',dFlag);
%
%   [mov,vObj] = ieMovie(rand(50,50,3,100),'step',2,'dFlag',dFlag);
%
% ISETBIO Team (BW) 2016


%% Parse inputs
p = inputParser;
p.addRequired('data',@isnumeric);
p.addParameter('vname','',@ischar);
p.addParameter('FrameRate',20,@isnumeric);
p.addParameter('step',1,@isnumeric);
p.addParameter('show',true,@islogical);
p.addParameter('gamma',1,@isnumeric);
p.KeepUnmatched = true;

p.parse(data,varargin{:});
data  = p.Results.data;
step  = p.Results.step;
show  = p.Results.show;
vname      = p.Results.vname;
FrameRate  = p.Results.FrameRate;
gamma      = p.Results.gamma;

%% Create the movie and video object

vObj = [];

% Could be monochrome or rgb
tDim = ndims(data);
nFrames = size(data, tDim);

% Scale and gamma correct mov
data = ieScale(data,0,1) .^ gamma;
mind = min(data(:)); maxd = max(data(:));
% A name for writing was passed
% So write and show the movie and write to file
if ~isempty(vname)
    axis image
    
    % When dFlag is a struct, show the move and save it in a file
    vObj = VideoWriter(vname);
    vObj.FrameRate = FrameRate;
    open(vObj);
    if isequal(tDim,4)
        % RGB data
        for ii = 1:step:nFrames
            imagesc(data(:,:,:,ii));
            caxis([mind maxd]);
            F = getframe; 
            writeVideo(vObj,F);
        end
    elseif isequal(tDim,3)
        % Monochrome data
        colormap(gray)
        for ii = 1:step:nFrames
            imagesc(data(:,:,ii));
            caxis([mind maxd]);
            F = getframe;
            writeVideo(vObj,F);
        end
    end
    close(vObj);

elseif show
    % Show, but don't write
    if isequal(tDim,4)
        % RGB data
        for ii=1:size(data,tDim)
            imagesc(data(:,:,:,ii)); axis image; caxis([mind maxd]); drawnow; 
        end
    elseif isequal(tDim,3)
        colormap(gray);
        for ii = 1:nFrames
            imagesc(data(:,:,ii)); axis image; caxis([mind maxd]); drawnow; 
        end
    end
end

end




