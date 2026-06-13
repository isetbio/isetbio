function mp4File = movie(obj,timeAxis, signal,varargin)
% Create a movie from the excitations or the photocurrent
%
% Synopsis
%   mp4File = cMosaic.movie(timeAxis, signal, varargin)
%
% Brief description
%   After computing a time series of cone excitations or photocurrent, send
%   in the signal and the time axis to create an mp4 movie.
%
% Inputs
%   obj - @cMosaic
%
% Optional key/value pairs
%   'file name'  -   filename with mp4 extension, default is 'test.mp4'
%   'photocurrent' - boolean, default is false
%
% Returns
%    file name
%

%% Inputs

varargin = ieParamFormat(varargin);
p = inputParser;
p.addRequired('obj',@(x)(isa(x,'cMosaic')));
p.addRequired('timeAxis',@isvector);
p.addRequired('signal',@(x)(ndims(x)==3));

p.addParameter('filename',fullfile(isetRootPath,'local','tmp.mp4'),@ischar);
p.addParameter('photocurrent',false);
p.addParameter('sample',1,@isinteger);

p.parse(obj,timeAxis,signal,varargin{:})

mp4File = p.Results.filename;
% if isempty(mp4File)
%     mp4File = fullfile(isetRootPath,'local','tmp.mp4');
% end

%%  Create the mp4 file

vidfile = VideoWriter(mp4File,'MPEG-4');
open(vidfile);

% You can play this video with VLC
for ii=1:numel(timeAxis)
    % Maybe this should be plot 'photocurrent' ....
    if ii == 1
        [~,hdl] = obj.plot('excitations',signal(p.Results.sample,ii,:));
    else
        obj.plot('excitations',signal(p.Results.sample,ii,:),'hdl',hdl);
    end
    drawnow;

    % Force the figure to be a fixed size so writeVideo never breaks
    if ii == 1
        set(hdl, 'Position', [100 100 800 600]);
        drawnow;
        frameSize = size(getframe(gcf).cdata);
    end
    thisImg = getframe(gcf);
    if sum(abs(size(thisImg.cdata) - frameSize)) > 0
        thisImg.cdata = imresize(thisImg.cdata, frameSize(1:2));
    end
    writeVideo(vidfile,thisImg);
end

close(vidfile)

end

