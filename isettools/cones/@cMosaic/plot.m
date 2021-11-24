function [uData, hdl] = plot(cmosaic,plotType, allE, varargin)
% plot methods for the cMosaic
%
% Syopsis
%    [uData, hdl] = plot(cmosaic, plotType, allE, varargin)
%
% Inputs
%    cmosaic   - cMosaic class
%    plotType  - See below.  Many.
%
% Optional key/val pairs
%    roi 
%    cone type
%
% Output
%    uData - Struct with the plot data
%    hdl   - Plot figure handle.  Use get(hdl,'CurrentAxes') for axis
%
% Plot types
%
%    excitations - pull out excitations of various types
%       horizontal line
%       vertical line (TBD)
%       roi - Does something for a line, rect or ellipse
%
%    roi - show the ROI superimposed on the axis
%
%
%
% See also
%   t_cMosaicBasic
%

%% Input parser

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cmosaic',@(x)(isa(x,'cMosaic')));
p.addRequired('plotType',@ischar);
p.addRequired('allE',@isnumeric);

% Excitations if precomputed
p.addParameter('conetype',{'l','m','s'},@(x)(ischar(x) || iscell(x)));
p.addParameter('roi',[],@(x)(isa(x,'regionOfInterest')));

% Horizontal line key val pairs
p.addParameter('ydeg',0,@isnumeric);
p.addParameter('thickness',0.1,@isnumeric);

p.parse(cmosaic,plotType,allE,varargin{:});

% Force cone type to a cell array
conetype = p.Results.conetype;
if ischar(conetype), conetype = {conetype}; 
end

%% Different types of plots

switch ieParamFormat(plotType)
    case {'excitations','activations'}
        % Show the activations in an image
        % We should choose one - excitations or activations - for
        % consistency. 

        params = cmosaic.visualize('params');
        params.activation = allE;
        
        % Return
        tmp = cmosaic.visualize(params);
        hdl = tmp.figureHandle;
        uData.allE = allE;
        
    case 'excitationshorizontalline'
        % Plot the excitations along an entire horizontal line in the
        % excitations
        %
        % You can select a single cone type.  Otherwise all three
        % are plotted in RGB colors.
        xminDeg = cmosaic.eccentricityDegs(1) - cmosaic.sizeDegs(1)/2;
        xmaxDeg = cmosaic.eccentricityDegs(1) + cmosaic.sizeDegs(1)/2;
        yDeg = p.Results.ydeg;
        
        if isempty(p.Results.roi)
            roi = regionOfInterest('shape', 'line', ...
                'from', [xminDeg,yDeg], 'to', [xmaxDeg,yDeg], ...
                'thickness', p.Results.thickness);
        else
            roi = p.Results.roi;
        end
        
        hdl = ieNewGraphWin;        
        for ii = 1:numel(conetype)
            [roiE, roiIdx] = cmosaic.excitations('roi',roi,...
                'conetype',conetype{ii},...
                'all excitations',allE);
            
            % The positions of the cones in the ROI
            pos = cmosaic.coneRFpositionsDegs(roiIdx,:);
            hold on;
            plot(pos(:,1),squeeze(roiE),[coneColor(conetype{ii}),'o']);
        end
        
        hold off; grid on
        str = sprintf('Excitations (%.1f ms)',cmosaic.integrationTime*1e3);
        xlabel('Horizontal position (deg)'); ylabel(str); 
        
        uData.pos = pos;
        uData.roiE = roiE;
        uData.roiIdx = roiIdx;
        
    case {'roi'}
        % Draw the ROI superimposed on the excitations
        fillColor = [0.2 0.3 1];

        % Excitations
        [~,hdl] = cmosaic.plot('excitations',allE);
        axesHandle = get(hdl,'CurrentAxes');
        
        roi = p.Results.roi;
        % Show the ROI as a horizontal line on the activations
        roiOutline = roi.outline();
        
        % Draw it but keeping the activations image
        hold(axesHandle, 'on');
        
        patch(axesHandle, roiOutline.x, roiOutline.y, ...
            fillColor, 'FaceAlpha', 0.5, 'EdgeColor', fillColor*0.5, 'LineWidth', 1.0);                        
        
        uData.excitations = allE;
        uData.roi = roi;
    case {'excitationsroi'}
        % Make a plot of the excitations in an ROI
        roi = p.Results.roi;
        switch roi.shape
            case 'rect'
                disp('Rect NYI')
            case 'ellipse'
                % Histogram of L,M,S???
                % 3d plot of L,M,S ???
                disp('ellipse NYI')
            case 'line'
                hdl = ieNewGraphWin;

                % Plot excitations along the line measured from 'from'                
                for ii = 1:numel(conetype)
                    [roiE, roiIdx] = cmosaic.excitations('roi',roi,...
                        'conetype',conetype{ii},...
                        'all excitations',allE);
                    pos = cmosaic.coneRFpositionsDegs(roiIdx,:);
                    linePosition = pos - roi.from;
                    lineDistance = sqrt(linePosition(:,1).^2 + linePosition(:,2).^2);
                    
                    % Set the color and collapse data w.r.t. the x-axis
                    hold on;
                    plot(lineDistance,squeeze(roiE),[coneColor(conetype{ii}),'o']);
                    
                    uData.lineDistance{ii} = lineDistance;
                    uData.roiE{ii} = roiE;
                end
                xlabel('Line position (deg - from)');
                ylabel(sprintf('Excitations %.1f',cmosaic.integrationTime));
                grid on;

            otherwise
                error('Unknown roi shape %s\n',roi.shape);
        end
        
    otherwise
        error('Unknown plot type %s\n',plotType);
end


end

function str = coneColor(conetype)
% Make this a routine
switch lower(conetype)
    case 'l'
        str = 'r';
    case 'm'
        str = 'g';
    case 's'
        str = 'b';
end
end
