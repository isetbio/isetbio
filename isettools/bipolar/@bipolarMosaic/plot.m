function hdl = plot(obj, pType, varargin)
% Plot data from a bipolar mosaic object
%
% Syntax:
%
%   hdl = bp.plot(plotType, varargin)
%
% Description:
%    Plot data from a bipolar mosaic object
%
% Input:
%    pType  - plot type
%
% Optional Key/Value Pairs:
%    gamma  - controls image display
%    pos    - positions to plot for time series
%
% Output:
%    hdl    - the bipolar mosaic plot according to the provided parameters
%
% Notes
%    * Maybe this should be obj.fig??? - question from below
%    * Need to address the other notes throughout the function
%    * Note: BW fix response center plot.

%% History:
% 5/2016 JRG, BW (c) isetbio team
%
%    10/18/17  jnm  Comments & Formatting
%    12/21/17  baw  Changed example.  Tested 'help'.  Added otherwise
%    case.

%% Examples:
%{
   % ETTBSkip - Example is broken. Remove this line when fixed.
   s_initRetina;
   bpMosaic = bpL.mosaic{1};
   bpMosaic.plot('spatial rf')
   bpMosaic.plot('mosaic');
   bpMosaic.plot('response center');
   bpMosaic.plot('response time series', 'pos', [5 5]);
   bpMosaic.plot('response image', 'gamma', 0.3);
   bpMosaic.plot('response movie');
%}

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;
p.FunctionName  = mfilename;
p.KeepUnmatched = true;
%%%
% Make key properties that can be set required arguments, and require
% values along with key names.
allowPlots = {...
    'help', ...
    'responsetimeseries', 'responsecenter', 'responsesurround', ...
    'responseimage', 'responsemovie', ...
    'spatialrf', 'surroundrf', 'mosaic'};
p.addRequired('pType', @(x) any(validatestring(ieParamFormat(x), ...
    allowPlots)));

p.addParameter('gamma', 1, @isscalar);
p.addParameter('pos', [], @ismatrix);

% Additional parameters are read in the case statements, below.
p.parse(pType, varargin{:});

%% Set the window.

if ~strcmpi(pType, 'help'), hdl = gcf; end
%vcNewGraphWin([], 'upperLeftBig');
sz = size(obj.responseCenter);

% Options
switch ieParamFormat(pType)
    case 'help'
        %% Case: Help
        fprintf('\nKnown %s plot types\n--------------\n', class(obj));
        for ii=2:length(allowPlots)
            fprintf('\t%s\n', allowPlots{ii});
        end
        hdl = [];
        return;
    case 'spatialrf'
        %% Case: Spatial RF
        %   @bipolarMosaic.plot('spatial rf')
        srf = obj.sRFcenter - obj.sRFsurround;
        sz = size(srf);
        if isequal(sz, [1, 1])
            disp('spatial rf is an impulse');
            return;
        end

        x = (1:sz(2)) - mean(1:sz(2));
        y = (1:sz(1)) - mean(1:sz(1));
        surf(x, y, srf); colormap(parula);
        xlabel('Cone samples'); zlabel('Responsivity')

    case 'surroundrf'
        %% Case: Surround RF
        %   @bipolarMosaic.plot('surround rf')
        srf = obj.sRFsurround;
        sz = size(srf);
        if isequal(sz, [1, 1])
            disp('spatial rf is an impulse');
            return;
        end

        x = (1:sz(2)) - mean(1:sz(2));
        y = (1:sz(1)) - mean(1:sz(1));
        surf(x, y, srf);
        colormap(parula);
        xlabel('Cone samples');
        zlabel('Responsivity')

    case {'mosaic'}
        %% Case: Mosaic
        %   @bipolarMosaic.plot('mosaic')
        %
        % Shows the RF Array, gets contour lines for mosaic RFs
        % The cell locations are specified with respect to the cone mosaic
        % input layer.  We would like to present them in terms of microns
        % on the cone mosaic surface.  So, we transform the cell locations
        % to microns.

        %%%
        % These are sampled w.r.t. the input mosaic.  We convert to microns
        center = obj.cellLocation;
        %%%
        % List the (x, y) positions on the grid and count how many
        center = reshape(center, [size(obj.cellLocation, 1) ...
            * size(obj.cellLocation, 2), 2]);
        nCells = size(center, 1);
        [r, c, ~] = size(obj.cellLocation);
        %%%
        % Convert sample grid positions to distance in microns
        metersPerBipolar = obj.patchSize ./ [r, c];
        %%%
        % Determine the RF radius.
        % Calculate the radius.  This seems too special case.  Ask JRG what
        % he intended here.  It seems like he thinks the radius is
        % predetermined to be 1.  But ...
        % At this moment, we are calculating the radius of the support in
        % microns. We probably want to have parameters that define the
        % support and the spread separately.
        center = 1e6 * center * diag(metersPerBipolar(:));
        % Centers in microns
        metersPerInput   = obj.input.patternSampleSize(1);
        radius = (1e6 * metersPerInput * size(obj.sRFcenter, 1)) / 2;
        %%%
        % At this point we should have centers and radius in terms of
        % microns.  If life is goo
        xMin = min(center(:, 2)) - 3 * radius;
        xMax = max(center(:, 2)) + 3 * radius;
        yMin = min(center(:, 2)) - 3 * radius;
        yMax = max(center(:, 2)) + 3 * radius;
        titleS = sprintf('RF positions and sizes');
        %%%
        % If there are a lot ... sub sample
        maxSamples = 500;
        if size(center, 1) > maxSamples
            center = center(randi(nCells, [maxSamples, 1]), :);
            titleS = sprintf(['Sampled RF positions and sizes '...
                '(%d of %d)'], maxSamples, size(obj.cellLocation(:), 1));
        end
        %%%
        % The whole ellipse thing isn't handled really, is it?  I mean, the
        % parameter is fixed here, and not set.
        ellipseMatrix = [1 1 0];
        ieShape('ellipse', ...
            'center', center, ...
            'radius', radius, ...
            'ellipseParameters', ellipseMatrix, ...
            'color', 'b');
        %%%
        % Sets the axis limits
        set(gca, 'xlim', [xMin, xMax], 'ylim', [yMin, yMax]);
        xlabel(sprintf('Distance (\\mum)'), 'fontsize', 14);
        ylabel(sprintf('Distance (\\mum)'), 'fontsize', 14);
        title(titleS);

    case{'responsecenter'}
        %% Case: Response Center
        %   @bipolarMosaic.plot('response center')
        %   @bipolarMosaic.plot('response center', 'pos', [1, 1]);
        %
        %
        % TODO - add cell location
        pos = p.Results.pos;
        if isempty(pos)
            %%%
            % Put position in the rows, time in the columns
            responseRS = RGB2XWFormat(obj.responseCenter);
        else
            nPos = size(pos, 1);
            responseRS = zeros(nPos, size(obj.responseCenter, 3));
            for ii=1:nPos
                responseRS(ii, :) = obj.responseCenter(pos(ii, 1), ...
                    pos(ii, 2), :);
            end
        end
        %%%
        % Show it
        vcNewGraphWin;
        plot(obj.timeStep * (1:sz(3)), responseRS');
        xlabel('Time (sec)'); ylabel('Response (AU)');
        title('Center response(s)'); grid on

    case{'responsesurround'}
        %% Case: Response Surround
        %   @bipolarMosaic.plot('response surround')
        %   responseRS = reshape(obj.responseSurround, sz(1)*sz(2), sz(3));
        pos = p.Results.pos;
        if isempty(pos)
            %%%
            % Put position in the rows, time in the columns
            responseRS = RGB2XWFormat(obj.responseSurround);
        else
            nPos = size(pos, 1);
            responseRS = zeros(nPos, size(obj.responseSurround, 3));
            for ii=1:nPos
                responseRS(ii, :) = obj.responseSurround(pos(ii, 1), ...
                    pos(ii, 2), :);
            end
        end

        vcNewGraphWin;
        plot(obj.timeStep * (1:sz(3)), responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar surround response'); grid on

    case{'responsetimeseries'}
        %% Case: Response Time Series
        %   @bipolarMosaic.plot('response time series', 'pos', ...)
        %
        % Open a new window and show the time series, but not in any
        % particular organized way. Only up to 200 samples are shown.
        % Random draws.  We should allow this to be controlled.
        %
        pos = p.Results.pos;
        if isempty(pos)
            responseRS = RGB2XWFormat(obj.responseCenter ...
                - obj.responseSurround);
            nCells = size(responseRS, 1);
        else
            nPos = size(pos, 1);
            responseRS = zeros(nPos, size(obj.responseSurround, 3));
            for ii=1:nPos
                responseRS(ii, :) = obj.responseCenter(pos(ii, 1), ...
                    pos(ii, 2), :) - obj.responseSurround(pos(ii, 1), ...
                    pos(ii, 2), :);
            end
        end
        tSamples = obj.timeStep * (1:size(obj.responseCenter, 3));

        vcNewGraphWin;
        if isempty(pos) &&  nCells > 100
            %%%
            % 100 randomly sampled positions
            cSamples = randi(nCells, [100, 1]);
            plot(tSamples, responseRS(cSamples, :));
            title('Bipolar current (100 samples)');
        elseif isempty(pos)
            %%%
            % All of the spatial samples
            plot(tSamples, responseRS);
            title('Bipolar current');
        else
            % From a single point
            plot(tSamples, responseRS)
            title(sprintf('Selected positions'));
        end

        xlabel(sprintf('Time (sec, \\Deltat %.0f ms)', ...
            obj.timeStep * 1e3));
        ylabel('Current (a.u.)');
        grid on

    case{'responseimage'}
        %% Case: Response Image
        %
        %   @bipolarMosaic.plot('response image')
        %
        response = (obj.get('response'));
        patchSizeUM = obj.patchSize * 1e6;
        dx = patchSizeUM / size(response, 1);  % Step in microns
        rowSamples = dx(1) * (1:size(response, 1));
        rowSamples = rowSamples - mean(rowSamples);
        colSamples = dx(2) * (1:size(response, 2));
        colSamples = colSamples - mean(colSamples);
        %%%
        % The gamma has to deal with negative numbers, sigh.
        if p.Results.gamma ~= 1
            img = mean(response, 3);
            img(img>=0) = img(img>=0).^p.Results.gamma;
            img(img<0)  = ((-1) * img(img<0) * p.Results.gamma) * -1;
        else
            img = mean(response, 3);
        end

        imagesc(rowSamples, colSamples, img);
        axis image; colormap(gray(256)); title('Current (a.u.)');
        xlabel(sprintf('Cell position (\\mum)'));

    case{'responsemovie'}
        %% Case: Response Movie
        % Pass the varargin along
        if ~isempty(varargin) && length(varargin) == 1
            %%%
            % Params are coded in a single struct
            varargin{1}.hf = hdl;
            if p.Results.gamma ~= 1
                varargin{1}.gamma = p.Results.gamma;
            end
            ieMovie(obj.get('response'), varargin{:});
        else
            %%%
            % Params are coded in param/value
            varargin{end+1} = 'gamma'; varargin{end+1} = p.Results.gamma;
            ieMovie(obj.get('response'), 'hf', hdl, varargin{:});
        end
    otherwise
        error('Unknown plot type %s\n',pType);
end

end
