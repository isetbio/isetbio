function [uData, hdl] = plot(cmosaic, plotType, allE, varargin)
% Plot common cMosaic views using a simple interface.
%
% Syntax
%   cmosaic.plot('mosaic')
%   cmosaic.plot('excitations', excitations)
%   cmosaic.plot('eye movement path')
%   cmosaic.plot('excitations and eye movements', excitations)
%
% Inputs
%   cmosaic  - A cMosaic object.
%   plotType - Named view to generate.
%   allE     - Cone excitations. May be a cone vector or a
%              trials-by-time-points-by-cones response array.
%
% Optional key/value pairs
%   trial         - Trial to display. Default is 1.
%   time point    - Time point to display, or 'last'. Default is 1.
%   roi           - regionOfInterest used by ROI views.
%   cone type     - Cone type or cell array of cone types.
%   xdeg, ydeg    - Position of vertical or horizontal response profiles.
%   thickness     - Profile ROI thickness in degrees.
%   figure handle - Figure in which to draw.
%   axes handle   - Axes in which to draw.
%   plot title    - Plot title.
%   label cones   - Label cone types in an excitation map.
%   data only     - Return plot data without drawing. Default is false.
%
% Outputs
%   uData - Struct containing selected and plotted data, including the axes.
%   hdl   - Figure handle, or empty when 'data only' is true.
%
% Description
%   This method provides common scientific views with sensible defaults.
%   Use cMosaic.visualize directly for specialized rendering controls.
%
% See also
%   cMosaic.visualize, cMosaic.excitations

if ~exist('allE', 'var'), allE = []; end

varargin = ieParamFormat(varargin);
plotType = ieParamFormat(plotType);

p = inputParser;
p.addRequired('cmosaic', @(x)isa(x, 'cMosaic'));
p.addRequired('plotType', @ischar);
p.addRequired('allE', @isnumeric);
p.addParameter('conetype', {'l', 'm', 's'}, @(x)ischar(x) || iscell(x));
p.addParameter('roi', [], @(x)isempty(x) || isa(x, 'regionOfInterest'));
p.addParameter('plottitle', '', @(x)ischar(x) || islogical(x) || iscell(x));
p.addParameter('dataonly', false, @islogical);
p.addParameter('labelcones', false, @islogical);
p.addParameter('lens', [], @(x)isempty(x) || isa(x, 'Lens'));
p.addParameter('trial', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('timepoint', 1, @(x)(isnumeric(x) && isscalar(x) && x >= 1) || ...
    (ischar(x) && strcmpi(x, 'last')));
p.addParameter('xdeg', 0, @isnumeric);
p.addParameter('ydeg', 0, @isnumeric);
p.addParameter('thickness', 0.1, @isnumeric);
p.addParameter('figurehandle', [], @(x)isempty(x) || isgraphics(x, 'figure'));
p.addParameter('axeshandle', [], @(x)isempty(x) || isgraphics(x, 'axes'));
% Backward-compatible figure-handle alias.
p.addParameter('hdl', [], @(x)isempty(x) || isgraphics(x, 'figure'));

% Time series parameters
p.addParameter('timeaxis', [], @isnumeric);
p.addParameter('coneindices', [], @isnumeric);
p.parse(cmosaic, plotType, allE, varargin{:});

figureHandle = p.Results.figurehandle;
if isempty(figureHandle), figureHandle = p.Results.hdl; end
axesHandle = p.Results.axeshandle;
if ~isempty(axesHandle), figureHandle = ancestor(axesHandle, 'figure'); end

conetype = p.Results.conetype;
if ischar(conetype), conetype = {conetype}; end

uData = struct('plotType', plotType, 'figureHandle', figureHandle, ...
    'axesHandle', axesHandle);
hdl = figureHandle;

switch plotType
    case {'conemosaic', 'mosaic'}
        if ~p.Results.dataonly
            % Preserve the original workhorse path for basic mosaic plots.
            % This direct call predates the generalized plot interface and
            % has been exercised extensively across MATLAB releases.
            visualizationParams = cmosaic.visualize;
            uData.figureHandle = visualizationParams.figureHandle;
            uData.axesHandle = visualizationParams.axesHandle;
            hdl = visualizationParams.figureHandle;
        end

    case {'excitations', 'activations'}
        % This view has historically been used for large mosaics.  If MATLAB
        % Desktop exits natively here, suspect a graphics/resource regression
        % rather than invalid excitation data.  Native exits bypass try/catch;
        % record the mosaic size and close residual MATLAB processes before
        % reproducing the problem in a fresh session.
        [selectedE, trial, timePoint] = localSelectExcitations(...
            cmosaic, allE, p.Results.trial, p.Results.timepoint);
        uData.excitations = selectedE;
        uData.trial = trial;
        uData.timePoint = timePoint;
        if ~p.Results.dataonly
            % Use the established parameter-struct invocation.  Besides
            % retaining long-tested rendering behavior, this still permits
            % the newer plot interface to select a trial/time point and to
            % supply existing figure or axes handles.
            visualizationParams = cmosaic.visualize('params');
            visualizationParams.activation = selectedE;
            visualizationParams.plotTitle = p.Results.plottitle;
            visualizationParams.verticalActivationColorBar = true;
            visualizationParams.figureHandle = figureHandle;
            visualizationParams.axesHandle = axesHandle;
            visualizationParams.labelConesInActivationMap = p.Results.labelcones;
            visualizationParams = cmosaic.visualize(visualizationParams);
            uData.figureHandle = visualizationParams.figureHandle;
            uData.axesHandle = visualizationParams.axesHandle;
            hdl = visualizationParams.figureHandle;
        end

    case {'eyemovementpath'}
        [trial, timePoints, currentEMposition] = localEyeMovementSelection(...
            cmosaic, p.Results.trial, p.Results.timepoint);
        uData.trial = trial;
        uData.timePoints = timePoints;
        uData.currentEMposition = currentEMposition;
        if ~p.Results.dataonly
            [uData, hdl] = localVisualize(cmosaic, [], p, true, ...
                figureHandle, axesHandle, uData);
        end

    case {'excitationsandeyemovements'}
        [selectedE, trial, timePoint] = localSelectExcitations(...
            cmosaic, allE, p.Results.trial, p.Results.timepoint);
        [~, timePoints, currentEMposition] = localEyeMovementSelection(...
            cmosaic, trial, timePoint);
        uData.excitations = selectedE;
        uData.trial = trial;
        uData.timePoint = timePoint;
        uData.timePoints = timePoints;
        uData.currentEMposition = currentEMposition;
        if ~p.Results.dataonly
            [uData, hdl] = localVisualize(cmosaic, selectedE, p, true, ...
                figureHandle, axesHandle, uData);
        end

    case {'timeseries', 'excitationstimeseries'}
        if isempty(allE)
            error('cMosaic:MissingExcitations', 'This plot type requires cone excitations.');
        end
        if isempty(p.Results.coneindices)
            if numel(size(allE)) == 3
                [~, maxIdx] = max(mean(allE, [1 2]));
            else
                [~, maxIdx] = max(mean(allE, 1));
            end
            [~,~,coneIndices] = ind2sub(size(allE), maxIdx);
            coneIndices = coneIndices(:)';
        else
            coneIndices = p.Results.coneindices(:)';
        end
        
        if ismatrix(allE)
            allE = reshape(allE, [1, size(allE,1), size(allE,2)]);
        end

        if isempty(p.Results.timeaxis)
            timeAxis = 1:size(allE, 2);
        else
            timeAxis = p.Results.timeaxis(:)';
        end

        uData.timeAxis = timeAxis;
        uData.excitations = allE(:, :, coneIndices);
        uData.coneIndices = coneIndices;

        if ~p.Results.dataonly
            [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle);
            hold(axesHandle, 'on');
            
            for cc = 1:length(coneIndices)
                cIdx = coneIndices(cc);
                if size(allE,1) > 1
                    timeSeries = squeeze(allE(:,:,cIdx));
                    if isvector(timeSeries), timeSeries = timeSeries(:)'; end
                    plot(axesHandle, timeAxis, timeSeries, 'b--', 'LineWidth', 0.5);
                end
                
                meanSeries = squeeze(mean(allE(:,:,cIdx),1));
                if isvector(meanSeries), meanSeries = meanSeries(:)'; end
                plot(axesHandle, timeAxis, meanSeries, 'k-', 'LineWidth', 2.0);
            end

            xlabel(axesHandle, 'Time');
            ylabel(axesHandle, 'Excitations');
            if ~isempty(p.Results.plottitle) && ischar(p.Results.plottitle)
                title(axesHandle, p.Results.plottitle);
            else
                title(axesHandle, 'Excitation Time Series');
            end
            grid(axesHandle, 'on');
            hold(axesHandle, 'off');
            hdl = figureHandle;
        end


    case {'excitationshorizontalline', 'excitationsverticalline'}
        [selectedE, trial, timePoint] = localSelectExcitations(...
            cmosaic, allE, p.Results.trial, p.Results.timepoint);
        if strcmp(plotType, 'excitationshorizontalline')
            limits = cmosaic.eccentricityDegs(1) + cmosaic.sizeDegs(1)/2*[-1 1];
            roi = localProfileROI(p.Results.roi, [limits(1) p.Results.ydeg], ...
                [limits(2) p.Results.ydeg], p.Results.thickness);
            positionDimension = 1;
            axisLabel = 'Horizontal position (deg)';
        else
            limits = cmosaic.eccentricityDegs(2) + cmosaic.sizeDegs(2)/2*[-1 1];
            roi = localProfileROI(p.Results.roi, [p.Results.xdeg limits(1)], ...
                [p.Results.xdeg limits(2)], p.Results.thickness);
            positionDimension = 2;
            axisLabel = 'Vertical position (deg)';
        end

        [uData, figureHandle, axesHandle] = localProfileData(cmosaic, ...
            selectedE, roi, conetype, positionDimension, p.Results.dataonly, ...
            figureHandle, axesHandle);
        uData.plotType = plotType;
        uData.trial = trial;
        uData.timePoint = timePoint;
        if ~p.Results.dataonly
            xlabel(axesHandle, axisLabel);
            ylabel(axesHandle, sprintf('Excitations (%.1f ms)', cmosaic.integrationTime*1e3));
            if ~isempty(p.Results.plottitle), title(axesHandle, p.Results.plottitle); end
        end
        hdl = figureHandle;

    case {'roi'}
        localRequireROI(p.Results.roi);
        [selectedE, trial, timePoint] = localSelectExcitations(...
            cmosaic, allE, p.Results.trial, p.Results.timepoint);
        uData.excitations = selectedE;
        uData.roi = p.Results.roi;
        uData.trial = trial;
        uData.timePoint = timePoint;
        if ~p.Results.dataonly
            [uData, hdl] = localVisualize(cmosaic, selectedE, p, false, ...
                figureHandle, axesHandle, uData);
            roiOutline = p.Results.roi.outline();
            hold(uData.axesHandle, 'on');
            patch(uData.axesHandle, roiOutline.x, roiOutline.y, [0.2 0.3 1], ...
                'FaceAlpha', 0.5, 'EdgeColor', [0.1 0.15 0.5], 'LineWidth', 1);
            hold(uData.axesHandle, 'off');
        end

    case {'excitationsroi'}
        localRequireROI(p.Results.roi);
        [selectedE, trial, timePoint] = localSelectExcitations(...
            cmosaic, allE, p.Results.trial, p.Results.timepoint);
        [uData, figureHandle, ~] = localROIData(cmosaic, selectedE, ...
            p.Results.roi, conetype, p.Results.dataonly, figureHandle, axesHandle);
        uData.plotType = plotType;
        uData.trial = trial;
        uData.timePoint = timePoint;
        hdl = figureHandle;

    case {'spectralqe'}
        if isempty(p.Results.lens)
            thisLens = Lens('wave', cmosaic.wave);
        else
            thisLens = p.Results.lens;
            thisLens.wave = cmosaic.wave;
        end
        uData.spectralQE = diag(thisLens.transmittance)*cmosaic.qe;
        if ~p.Results.dataonly
            [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle);
            plot(axesHandle, cmosaic.wave, uData.spectralQE, 'LineWidth', 2);
            xlabel(axesHandle, 'Wavelength (nm)');
            ylabel(axesHandle, 'Spectral quantum efficiency');
            grid(axesHandle, 'on');
            if ~isempty(p.Results.plottitle), title(axesHandle, p.Results.plottitle); end
        end
        uData.figureHandle = figureHandle;
        uData.axesHandle = axesHandle;
        hdl = figureHandle;

    case {'pigmentquantalefficiency'}
        uData.quantalEfficiency = cmosaic.pigment.quantalEfficiency;
        if ~p.Results.dataonly
            [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle);
            plot(axesHandle, cmosaic.wave, uData.quantalEfficiency, 'LineWidth', 2);
            xlabel(axesHandle, 'Wavelength (nm)');
            ylabel(axesHandle, 'Quantal efficiency');
            if ~isempty(p.Results.plottitle), title(axesHandle, p.Results.plottitle); end
        end
        uData.figureHandle = figureHandle;
        uData.axesHandle = axesHandle;
        hdl = figureHandle;

    otherwise
        error('cMosaic:UnknownPlotType', 'Unknown plot type ''%s''.', plotType);
end

if p.Results.dataonly, hdl = []; end

end

function [uData, figureHandle] = localVisualize(cmosaic, activation, p, ...
    showEyeMovements, figureHandle, axesHandle, uData)
args = {'figureHandle', figureHandle, 'axesHandle', axesHandle, ...
    'plotTitle', p.Results.plottitle};
if ~isempty(activation)
    args = [args {'activation', activation, 'verticalActivationColorBar', true, ...
        'labelConesInActivationMap', p.Results.labelcones}];
end
if showEyeMovements
    args = [args {'currentEMposition', uData.currentEMposition, ...
        'displayedEyeMovementData', struct('trial', uData.trial, ...
        'timePoints', uData.timePoints)}];
end
params = cmosaic.visualize(args{:});
figureHandle = params.figureHandle;
uData.figureHandle = params.figureHandle;
uData.axesHandle = params.axesHandle;
end

function [selectedE, trial, timePoint] = localSelectExcitations(cmosaic, allE, trial, timePoint)
if isempty(allE)
    error('cMosaic:MissingExcitations', 'This plot type requires cone excitations.');
end
if isvector(allE)
    selectedE = allE(:);
    trial = 1;
    timePoint = 1;
else
    if size(allE, 3) ~= cmosaic.conesNum
        error('cMosaic:ExcitationSizeMismatch', ...
            'The third excitation dimension must equal the number of cones (%d).', cmosaic.conesNum);
    end
    if ischar(timePoint), timePoint = size(allE, 2); end
    if trial > size(allE, 1) || timePoint > size(allE, 2)
        error('cMosaic:ExcitationIndexOutOfRange', ...
            'Requested trial %d and time point %d exceed the excitation dimensions.', trial, timePoint);
    end
    selectedE = squeeze(allE(trial, timePoint, :));
end
if numel(selectedE) ~= cmosaic.conesNum
    error('cMosaic:ExcitationSizeMismatch', ...
        'Excitations must contain one value for each of the %d cones.', cmosaic.conesNum);
end
end

function [trial, timePoints, currentEMposition] = localEyeMovementSelection(cmosaic, trial, timePoint)
if isempty(cmosaic.fixEMobj)
    error('cMosaic:MissingEyeMovements', ...
        'Generate eye movements with cMosaic.emGenSequence before plotting them.');
end
nTrials = size(cmosaic.fixEMobj.emPosArcMin, 1);
nTimePoints = size(cmosaic.fixEMobj.emPosArcMin, 2);
if ischar(timePoint), timePoint = nTimePoints; end
if trial > nTrials || timePoint > nTimePoints
    error('cMosaic:EyeMovementIndexOutOfRange', ...
        'Requested trial %d and time point %d exceed the eye-movement dimensions.', trial, timePoint);
end
timePoints = 1:timePoint;
currentEMposition = squeeze(cmosaic.fixEMobj.emPosArcMin(trial, timePoint, :))/60;
end

function roi = localProfileROI(roi, from, to, thickness)
if isempty(roi)
    roi = regionOfInterest('shape', 'line', 'from', from, 'to', to, ...
        'thickness', thickness);
end
end

function localRequireROI(roi)
if isempty(roi)
    error('cMosaic:MissingROI', 'This plot type requires a regionOfInterest.');
end
end

function [uData, figureHandle, axesHandle] = localProfileData(cmosaic, allE, roi, ...
    conetype, positionDimension, dataOnly, figureHandle, axesHandle)
uData = struct('roi', roi, 'positions', {cell(numel(conetype), 1)}, ...
    'roiExcitations', {cell(numel(conetype), 1)}, 'roiIndices', {cell(numel(conetype), 1)});
if ~dataOnly, [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle); end
for ii = 1:numel(conetype)
    [roiE, roiIdx] = cmosaic.excitations('roi', roi, ...
        'conetype', conetype{ii}, 'all excitations', allE);
    positions = cmosaic.coneRFpositionsDegs(roiIdx, :);
    uData.positions{ii} = positions;
    uData.roiExcitations{ii} = roiE;
    uData.roiIndices{ii} = roiIdx;
    if ~dataOnly
        hold(axesHandle, 'on');
        thisPlot = plot(axesHandle, positions(:, positionDimension), squeeze(roiE), ...
            [localConeColor(conetype{ii}) 'o']);
        set(thisPlot, 'MarkerFaceColor', localConeColor(conetype{ii}));
    end
end
if ~dataOnly
    hold(axesHandle, 'off');
    grid(axesHandle, 'on');
end
uData.figureHandle = figureHandle;
uData.axesHandle = axesHandle;
end

function [uData, figureHandle, axesHandle] = localROIData(cmosaic, allE, roi, ...
    conetype, dataOnly, figureHandle, axesHandle)
uData = struct('roi', roi, 'roiExcitations', {cell(numel(conetype), 1)}, ...
    'roiIndices', {cell(numel(conetype), 1)});
if ~dataOnly, [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle); end
for ii = 1:numel(conetype)
    [roiE, roiIdx] = cmosaic.excitations('roi', roi, ...
        'conetype', conetype{ii}, 'all excitations', allE);
    uData.roiExcitations{ii} = roiE;
    uData.roiIndices{ii} = roiIdx;
    if strcmp(roi.shape, 'line')
        positions = cmosaic.coneRFpositionsDegs(roiIdx, :);
        linePosition = positions - roi.from;
        uData.lineDistance{ii} = sqrt(sum(linePosition.^2, 2));
    end
    if ~dataOnly
        hold(axesHandle, 'on');
        if strcmp(roi.shape, 'line')
            plot(axesHandle, uData.lineDistance{ii}, squeeze(roiE), ...
                [localConeColor(conetype{ii}) 'o']);
        else
            histogram(axesHandle, roiE, 'FaceColor', localConeColor(conetype{ii}), ...
                'EdgeColor', localConeColor(conetype{ii}), 'NumBins', 20);
        end
    end
end
if ~dataOnly
    hold(axesHandle, 'off');
    grid(axesHandle, 'on');
    if strcmp(roi.shape, 'line')
        xlabel(axesHandle, 'Line position (deg from start)');
    else
        xlabel(axesHandle, 'Excitations');
    end
    ylabel(axesHandle, 'Excitations');
end
uData.figureHandle = figureHandle;
uData.axesHandle = axesHandle;
end

function [figureHandle, axesHandle] = localAxes(figureHandle, axesHandle)
if isempty(axesHandle)
    if isempty(figureHandle)
        figureHandle = ieFigure;
        axesHandle = gca;
    else
        axesHandle = axes('Parent', figureHandle);
    end
else
    figureHandle = ancestor(axesHandle, 'figure');
    cla(axesHandle);
end
end

function color = localConeColor(conetype)
switch lower(conetype)
    case 'l'
        color = 'r';
    case 'm'
        color = 'g';
    case 's'
        color = 'b';
    otherwise
        error('cMosaic:UnknownConeType', 'Unknown cone type ''%s''.', conetype);
end
end
