function [uData, hf] = plot(obj, type, varargin)
% Plot function for coneMosaicHex (subclass of coneMosaic)
%
%    [uData, hf] = coneMosaicHex.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure handle or control structure, the meaning of value is%             []            - create plot in new figure
%             figure handle - plot in figure specified by hf
%             'none'        - don't plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
% See coneMosaic.plot for plot types
%
% Example:
%    rgc.mosaic{1}.plot(type)
%
% HJ/BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj');
p.addRequired('type', @isstr);               % Type of plot

p.addParameter('hf', []);                    % figure handle
p.addParameter('oi',[],@isstruct);           % Used for spectral qe

p.parse(obj,type, varargin{:});
hf = p.Results.hf;

uData = [];

%% Figure and axes
if isempty(hf), hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure'), figure(hf);
elseif isgraphics(hf, 'axes'), axes(hf);
end

% set color order so that LMS plots as RGB
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');
    if isgraphics(hf,'axes')
        set(get(hf,'parent'),'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
    else  % Figure
        set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
    end
end

%% Pick up the specialized hex routines here

% If not one of ours, pass on to the base class in the otherwise case.
switch ieParamFormat(type);
    case 'conemosaic'
        % Try the NP plotting routine from here ...
        % It brings the image up in the wrong window because we didn't pass
        % the hf and set the axis.  But the good news is, it gets here and
        % runs.  Next, on to fixing the varargin{} and such.
        plotHexMosaic(obj);  % Default arguments for now
    case 'meanabsorptions'

        % Copied from t_coneMosaicHex3.m, line 70
        isomerizationsHex = obj.absorptions;
        nonNullConeIndices = obj.pattern > 1;
        isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
                
        % Render activation images for the hex mosaic
        [activationsHexImage,activationsType,supX,supY] = obj.computeActivationDensityMap(isomerizationsHex);
                
        % Display results
        % imagesc(1:obj.cols, 1:obj.rows, mean(activationsHexImage,3));
        imagesc(supX,supY,mean(activationsHexImage,3));
        % set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
        axis 'image'; axis 'xy'
        colormap gray;
        
        hold on;
        sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
        sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);
        dx = obj.pigment.pdWidth;
%         hold off
        axis 'equal'; axis 'xy'
        xTicks = [sampledHexMosaicXaxis(1) obj.center(1) sampledHexMosaicXaxis(end)];
        yTicks = [sampledHexMosaicYaxis(1) obj.center(2) sampledHexMosaicYaxis(end)];
        xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
        yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
        set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], 'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
        box on; grid off;
        set(gca, 'XLim', [sampledHexMosaicXaxis(1)-dx sampledHexMosaicXaxis(end)+dx]);
        set(gca, 'YLim', [sampledHexMosaicYaxis(1)-dx sampledHexMosaicYaxis(end)+dx]);
        drawnow

    case 'movieabsorptions'

        % Copied from t_coneMosaicHex3.m, line 70
        isomerizationsHex = obj.absorptions;
        nonNullConeIndices = obj.pattern > 1;
        isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
                
        % Render activation images for the hex mosaic
        [activationsHexImage,activationsType,supX,supY] = obj.computeActivationDensityMap(isomerizationsHex);
                
        activationsHexMovie = zeros([size(activationsHexImage),size(isomerizationsHex,3)]);
        for frameIndex = 1:size(isomerizationsHex,3)
            [activationsHexImage, ~] = obj.computeActivationDensityMap(isomerizationsHex(:,:,frameIndex));
            activationsHexMovie(:,:,frameIndex) = activationsHexImage;
        end
        
        % set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
        axis 'image'; axis 'xy'
        % title('hex mosaic isomerizations (all cones)', 'FontSize', 16);
        % % % % % % %
        uData = ieMovie(activationsHexMovie,varargin{:});
        
    case 'meancurrent'
        disp('NYI');
%         if isempty(obj.current)
%             error('no photocurrent data computed');
%         end
%         uData = mean(obj.current, 3);
%         if ~isequal(hf, 'none')
%             imagesc(uData); axis off; colorbar;
%             % title('Mean photocurrent (pA)');
%         end
%         colormap(gray); % Shows a numerical value
%         axis image;


        % Copied from t_coneMosaicHex3.m, line 70
        clear isomerizationsHex
        isomerizationsHex = obj.current;
        nonNullConeIndices = obj.pattern > 1;
        isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
                
        % Render activation images for the hex mosaic
        [activationsHexImage, activationsType,supX,supY] = obj.computeActivationDensityMap(isomerizationsHex);
                
        % Display results
        % imagesc(1:obj.cols, 1:obj.rows, mean(activationsHexImage,3));
        imagesc(supX,supY,mean(activationsHexImage,3));
        % set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
        axis 'image'; axis 'xy'
        colormap gray;
        
        hold on;
        sampledHexMosaicXaxis = obj.patternSupport(1,:,1) + obj.center(1);
        sampledHexMosaicYaxis = obj.patternSupport(:,1,2) + obj.center(2);
        dx = obj.pigment.pdWidth;
%         hold off
        axis 'equal'; axis 'xy'
        xTicks = [sampledHexMosaicXaxis(1) obj.center(1) sampledHexMosaicXaxis(end)];
        yTicks = [sampledHexMosaicYaxis(1) obj.center(2) sampledHexMosaicYaxis(end)];
        xTickLabels = sprintf('%2.0f um\n', xTicks*1e6);
        yTickLabels = sprintf('%2.0f um\n', yTicks*1e6);
        set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(gca, 'FontSize', 16, 'XColor', [0.1 0.2 0.9], 'YColor', [0.1 0.2 0.9], 'LineWidth', 1.0);
        box on; grid off;
        set(gca, 'XLim', [sampledHexMosaicXaxis(1)-dx sampledHexMosaicXaxis(end)+dx]);
        set(gca, 'YLim', [sampledHexMosaicYaxis(1)-dx sampledHexMosaicYaxis(end)+dx]);
%         drawnow

    case 'moviecurrent'
        
        disp('NYI');
%         % Current movie in gray scale
%         if isempty(obj.current)
%             if isempty(p.Results.hf), close(hf); end
%             error('no current data');
%         end
%         % Additional arguments may be the video file name, step, and
%         % FrameRate
%         uData = ieMovie(obj.current,varargin{:});
        

        % Copied from t_coneMosaicHex3.m, line 70
        isomerizationsHex = obj.current;
        nonNullConeIndices = obj.pattern > 1;
        isomerizationsRange = prctile(isomerizationsHex(:), [5 95]);
                
        % Render activation images for the hex mosaic
        [activationsHexImage, ~] = obj.computeActivationDensityMap(isomerizationsHex);
        
        activationsHexMovie = zeros([size(activationsHexImage),size(isomerizationsHex,3)]);
        for frameIndex = 1:size(isomerizationsHex,3)
            [activationsHexImage, ~] = obj.computeActivationDensityMap(isomerizationsHex(:,:,frameIndex));
            activationsHexMovie(:,:,frameIndex) = activationsHexImage;
        end
        
        % set(gca, 'CLim', isomerizationsRange, 'XTick', [], 'YTick', []);
        axis 'image'; axis 'xy'
        % title('hex mosaic isomerizations (all cones)', 'FontSize', 16);
        % % % % % % %
        uData = ieMovie(activationsHexMovie,varargin{:});
    otherwise
        % Not one of the hex image types.  So, pass up to the base class
        [uData,hf] = plot@coneMosaic(obj,type,varargin{:});
end
