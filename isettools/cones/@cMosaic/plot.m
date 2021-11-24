function [uData, hdl] = plot(cmosaic,plotType, allE, varargin)
% plot methods for the cMosaic
%
% Syopsis
%    [uData, hdl] = plot(cmosaic,plotType, allE, varargin)
%
% Inputs
%    cmosaic
%    plotType  - Horizontal line, ellipse region ...
%
% Optional key/val pairs
%    Many, grouped by plot type
%
% Output
%    uData - Struct with the plot data
%
% See also
%

%% Input parser

varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cmosaic',@(x)(isa(x,'cMosaic')));
p.addRequired('plotType',@ischar);
p.addRequired('allE',@isnumeric);

% Excitations if precomputed
p.addParameter('allexcitations',[],@isnumeric);
p.addParameter('conetype','',@ischar);
p.addParameter('roi',[],@(x)(isa(x,'regionOfInterest')));

% Horizontal line key val pairs
p.addParameter('ydeg',0,@isnumeric);
p.addParameter('thickness',0.1,@isnumeric);

p.parse(cmosaic,plotType,allE,varargin{:});

%% Different types of plots

switch ieParamFormat(plotType)
    case 'horizontalline'
        ct = {'l','m','s'};
        cOrder = {'r','g','b'};
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

        % Extract the excitations in the ROI
        if ~isempty(p.Results.conetype)
            [roiE, roiIdx] = cmosaic.excitations('roi',roi,...
                'conetype',p.Results.conetype,...
                'all excitations',allE);
            pos = cmosaic.coneRFpositionsDegs(roiIdx,:);
            
            % Set the color and collapse data w.r.t. the x-axis          
            [~,idx] = ismember(p.Results.conetype,ct);
            plot(pos(:,1),squeeze(roiE),[cOrder{idx},'o']);
            
        else
            for ii = 1:3
                [roiE, roiIdx] = cmosaic.excitations('roi',roi,...
                    'conetype',ct{ii},...
                    'all excitations',allE);
                
                % The positions of the cones in the ROI
                pos = cmosaic.coneRFpositionsDegs(roiIdx,:);
                hold on;
                plot(pos(:,1),squeeze(roiE),[cOrder{ii},'o']);
                
            end
        end
        
        hold off; grid on
        xlabel('Horizontal position (deg)');
        ylabel('Excitations'); 
        uData.pos = pos;
        uData.roiE = roiE;
        uData.roiIdx = roiIdx;
        
    otherwise
        error('Unknown plot type %s\n',plotType);
end


end
