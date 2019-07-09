function [uData, hf] = plot(obj, type, varargin)
% Plot function for coneMosaicHex (subclass of coneMosaic)
%
% Syntax:
%	[uData, hf] = coneMosaicHex.plot(type, [varargin])
%   [uData, hf] = plot(obj, type, [varargin])
%
% Inputs:
%    obj  - The cone mosaic hex object
%    type - A string denoting the type of the plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
% Optional key/value pairs:
%   'hf'  - figure handle or control structure, the meaning of value is
%               []            - create plot in new figure
%               figure handle - plot in figure specified by hf
%               'none'        - don't plot
%
% See Also:
%    coneMosaic.plot for plot types
%

% History:
%    xx/xx/16  HJ/BW  ISETBIO TEAM, 2016
%    02/20/18  jnm    Formatting
%    04/07/18  dhb    Skip broken example.

% Examples:
%{
    % ETTBSkip. This needs a mosaic defined before it could possibly work.
    rgc.mosaic{1}.plot(type)
%}

%% parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj');
p.addRequired('type', @isstr);        % Type of plot

p.addParameter('hf', []);             % figure handle
p.addParameter('oi', [], @isstruct);  % Used for spectral qe

p.parse(obj, type, varargin{:});
hf = p.Results.hf;

uData = [];

%% Figure and axes
if isempty(hf)
    hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure')
    figure(hf);
elseif isgraphics(hf, 'axes')
    axes(hf);
end

% Set color order so that LMS plots as RGB
% Matlab default is 7 colors, and we reorder
% If the user has changed the default, we leave it alone.
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');  
    if size(co,1) == 7    % Figure        
        if isgraphics(hf, 'axes')
            set(get(hf, 'parent'), ...
                'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
        else
            set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
        end
    end
end

%% Pick up the specialized hex routines here
% If not one of ours, pass on to the base class in the otherwise case.
switch ieParamFormat(type)
    case 'conemosaic'
        % Try the NP plotting routine from here ...
        % It brings the image up in the wrong window because we didn't pass
        % the hf and set the axis. But the good news is, it gets here and
        % runs. Next, on to fixing the varargin{} and such.
        %
        % [Note: DHB - I'm not sure I had to change this, but passing
        % varargin{:} onwards is a bit confusing to me, as hf can get
        % passed to the next routine twice if you're not careful. Probably
        % should explicitly set anything that we want plotHexMosaic to
        % inherit here.]

        % Default arguments for now
        % plotHexMosaic(obj, 'hf', hf, varargin{:});
        plotHexMosaic(obj, 'hf', hf);

    case 'meanabsorptions'
        plotHexMeanImage(obj, 'type', 'absorptions');
%                 clear isomerizationsHex
%         isomerizationsHex = obj.current;
%         % nonNullConeIndices = obj.pattern > 1;
%         % isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
%
%         % Render activation images for the hex mosaic
%         [activationsHexImage, ~, supX, supY] = ...
%              obj.computeActivationDensityMap(...
%              mean(isomerizationsHex(:, :, end), 3));
%
%         % Display results
%         % imagesc(1:obj.cols, 1:obj.rows, mean(activationsHexImage, 3));
%         imagesc(supX, supY, activationsHexImage);
%         % set(gca, 'CLim', isomerizationsRange, 'XTick', [], ...
%         %      'YTick', []);
%         axis 'image';
%         axis 'xy'
%         colormap gray;
%
%         hold on;
%         sampledHexMosaicXaxis = obj.patternSupport(1, :, 1) + ...
%               obj.center(1);
%         sampledHexMosaicYaxis = obj.patternSupport(:, 1, 2) + ...
%               obj.center(2);
%         dx = obj.pigment.pdWidth;
%
%         axis 'equal';
%         axis 'xy'
%         xTicks = [sampledHexMosaicXaxis(1), obj.center(1), ...
%             sampledHexMosaicXaxis(end)];
%         yTicks = [sampledHexMosaicYaxis(1), obj.center(2), ...
%             sampledHexMosaicYaxis(end)];
%         xTickLabels = sprintf('%2.0f um\n', xTicks * 1e6);
%         yTickLabels = sprintf('%2.0f um\n', yTicks * 1e6);
%         set(gca, 'XTick', xTicks, 'YTick', yTicks, ...
%             'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
%         set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], ...
%             'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
%         box on;
%         grid off;
%         isomerizationsRange = prctile(activationsHexImage(:), [5 95]);
%
%         set(gca, 'CLim', isomerizationsRange);
%         set(gca, 'XLim', [sampledHexMosaicXaxis(1) - dx ...
%             sampledHexMosaicXaxis(end) + dx]);
%         set(gca, 'YLim', [sampledHexMosaicYaxis(1) - dx ...
%             sampledHexMosaicYaxis(end) + dx]);

    case 'movieabsorptions'
        uData = movieHex(obj, 'type', 'absorptions');
%         % Copied from t_coneMosaicHex3.m, line 70
%         isomerizationsHex = obj.absorptions;
%         % nonNullConeIndices = obj.pattern > 1;
%         % isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
%
%         % Render activation images for the hex mosaic
%         [activationsHexImage, ~, ~, ~] = ...
%             obj.computeActivationDensityMap(isomerizationsHex);
%
%         activationsHexMovie = zeros([size(activationsHexImage), ...
%             size(isomerizationsHex, 3)]);
%         for frameIndex = 1:size(isomerizationsHex, 3)
%             [activationsHexImage, ~] = ...
%                 obj.computeActivationDensityMap(...
%                 isomerizationsHex(:, :, frameIndex));
%             activationsHexMovie(:, :, frameIndex) = activationsHexImage;
%         end
%
%         axis 'image';
%         axis 'xy'
%         uData = ieMovie(activationsHexMovie, varargin{:});

    case 'meancurrent'
        plotHexMeanImage(obj, 'type', 'current');

    case 'moviecurrent'
        uData = movieHex(obj, 'type', 'current');
%         % Copied from t_coneMosaicHex3.m, line 70
%         isomerizationsHex = obj.current;
%         % nonNullConeIndices = obj.pattern > 1;
%         % isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
%
%         % Render activation images for the hex mosaic
%         [activationsHexImage, ~] = obj.computeActivationDensityMap(...
%             isomerizationsHex);
%
%         activationsHexMovie = zeros([size(activationsHexImage), ...
%             size(isomerizationsHex, 3)]);
%         for frameIndex = 1:size(isomerizationsHex, 3)
%             [activationsHexImage, ~] = ...
%                 obj.computeActivationDensityMap(...
%                 isomerizationsHex(:, :, frameIndex));
%             activationsHexMovie(:, :, frameIndex) = activationsHexImage;
%         end
%
%         % set(gca, 'CLim', isomerizationsRange, 'XTick', [], ...
%         %    'YTick', []);
%         axis 'image';
%         axis 'xy'
%         % title('hex mosaic isomerizations (all cones)', 'FontSize', 16);
%         % % % % % % %
%         uData = ieMovie(activationsHexMovie, varargin{:});

    otherwise
        % Not one of the hex image types. So, pass up to the base class
        [uData, hf] = plot@coneMosaic(obj, type, varargin{:});
end
