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

p.addParameter('filename', fullfile(isetRootPath,'local','tmp.mp4'), ...
    @(x)(ischar(x) || (isstring(x) && isscalar(x))));
p.addParameter('photocurrent', false, @islogical);
p.addParameter('sample', 1, ...
    @(x)(isnumeric(x) && isscalar(x) && isfinite(x) && ...
    (x >= 1) && (round(x) == x)));

p.parse(obj,timeAxis,signal,varargin{:})

mp4File = p.Results.filename;
% if isempty(mp4File)
%     mp4File = fullfile(isetRootPath,'local','tmp.mp4');
% end

%%  Create the mp4 file

vidfile = VideoWriter(mp4File, 'MPEG-4');
open(vidfile);
writerCleanup = onCleanup(@() closeVideoWriter(vidfile));

% You can play this video with VLC
hdl = [];
axesHandle = [];
figureCleanup = [];
frameSize = [];
for ii = 1:numel(timeAxis)
    % Maybe this should be plot 'photocurrent' ....
    if ii == 1
        [plotData, hdl] = obj.plot('excitations', ...
            signal(p.Results.sample,ii,:));
        axesHandle = plotData.axesHandle;

        % cMosaic plotting commonly creates an ieFigure whose units are
        % normalized. Switch explicitly to pixels before assigning a
        % pixel-sized position, then restore the incoming unit state.
        originalFigureUnits = hdl.Units;
        hdl.Units = 'pixels';
        set(hdl, ...
            'Position', [100 100 800 600], ...
            'Resize', 'off');
        hdl.Units = originalFigureUnits;
        figureCleanup = onCleanup(@() closeVideoFigure(hdl));
    else
        % Reuse both handles. Passing only the figure causes visualize to
        % clear the figure, create new axes, and reset the figure geometry.
        obj.plot('excitations', signal(p.Results.sample,ii,:), ...
            'figureHandle', hdl, ...
            'axesHandle', axesHandle);
    end

    % Capture the returned figure, not whichever figure happens to be
    % current. getframe performs the required graphics update.
    thisImg = getframe(hdl);
    if isempty(frameSize)
        frameSize = size(thisImg.cdata);
    elseif ~isequal(size(thisImg.cdata), frameSize)
        % Resize defensively if the window manager changes drawable pixels.
        thisImg.cdata = imresize(thisImg.cdata, frameSize(1:2));
    end
    writeVideo(vidfile, thisImg);
end

close(vidfile);
clear writerCleanup;
clear figureCleanup;

end

function closeVideoWriter(vidfile)
% Close a partially written video after a MATLAB-level error.

try
    close(vidfile);
catch
    % The writer may already be closed.
end
end

function closeVideoFigure(hFig)
% Release the graphics resources used to render video frames.

if isgraphics(hFig)
    close(hFig);
end
end
