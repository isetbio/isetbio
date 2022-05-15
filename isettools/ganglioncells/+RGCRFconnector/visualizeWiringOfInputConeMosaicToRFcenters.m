function [hFig, ax, XLims, YLims] = visualizeWiringOfInputConeMosaicToRFcenters(...
    RGCRFinputs, RGCRFweights, allConePositions, allConeSpacings, allConeTypes, varargin)
% Visualize the wiring of the inputConeMosaic to the RGC RF centers
%
% Syntax:
%   [hFig, ax, XLims, YLims] =
%   RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
%      RGCRFinputs, RGCRFweights, allConePositions, allConeSpacings, ...
%      allConeTypes, varargin);
%
% Description:
%   Visualize the wiring of the inputConeMosaic to the RGC RF centers
%
% Inputs:
%    RGCRFinputs        - Cell array with indices of the input cones, one cell per each target RGC RF center
%    RGCRFweights       - Cell array with weights of the input cones, one cell per each target RGC RF center 
%    allConePositions   - [N x 2] matrix of (x,y) positions of the N cones in the mosaic
%    allConeSpacings    - [N x 1] vector of spacings of the N cones in the mosaic
%    allConeTypes       - [N x 1] vector of types of the N cones in the mosaic
%
% Outputs:
%    hFig               - The figure handle
%    ax                 - The axes handle
%    XLims, YLims       - The limits [min, max] of the x- and y-axes
%
% Optional key/value pairs
%   'figureHandle'                  - The figure handle on which to render the figure
%   'axesHandle'                    - The axes handle on which to render the figure
%   'XLims', 'YLims'                - The limits [min max] of the x- and y-axes
%   'thetaSamples'                  - How many samples to use for rendering the RF disk
%   'superimposeConeInputWiring'    - Logical. Whether to superimpose lines
%                                       from input cone to RGCRF centroid
%   'displayRGCID'                  - Logical. Whether to display the RGCID
%   'titleString'                   - A title string
%   'pauseAfterEachRGCisRendered'   - Pause for debug reasons
%   'videoOBJ'                      - VideoOBJ for exporting each connection step as a frame in a video
%
% History:
%   5/11/2022       NPC     Wrote it
%

    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('thetaSamples', 15, @isnumeric);
    p.addParameter('superimposeConeInputWiring', false, @islogical);
    p.addParameter('titleString', '', @(x)(isempty(x) || (ischar(x))));
    p.addParameter('pauseAfterEachRGCisRendered', false, @islogical);
    p.addParameter('displayRGCID', false, @islogical);
    p.addParameter('videoOBJ', [], @(x)(isempty(x)||isa(x, 'VideoWriter')));
    p.addParameter('XLims', [], @isnumeric);
    p.addParameter('YLims', [], @isnumeric);
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    titleString = p.Results.titleString;
    thetaSamples = p.Results.thetaSamples;
    superimposeConeInputWiring = p.Results.superimposeConeInputWiring;
    pauseAfterEachRGCisRendered = p.Results.pauseAfterEachRGCisRendered;
    displayRGCID = p.Results.displayRGCID;
    videoOBJ = p.Results.videoOBJ;
   
    XLims = p.Results.XLims;
    YLims = p.Results.YLims;

    % Generate disk outline
    thetas = linspace(0,360,thetaSamples);
    diskOutline = 0.5*[cosd(thetas); sind(thetas)]';


    % Initialize figure
    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1]);
        end
        ax = subplot('Position', [0.1 0.1 0.8 0.8]);
    end

    % Render cones
    hold(ax, 'on');
    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    coneColors = [1 0 0.2; ...
                  0 0.8 0.2; ...
                  0 0.5 1];
        
    for iConeType = 1:numel(coneTypes)
        idx = find(allConeTypes == coneTypes(iConeType));
        [f,v] = facesAndVertices(allConePositions(idx,:), 1.5*0.204*sqrt(2.0)*allConeSpacings(idx), diskOutline);
        theColor = squeeze(coneColors(iConeType,:));
        patch(ax,'Faces', f, 'Vertices', v, 'LineWidth', 0.5, 'FaceColor', theColor, 'EdgeColor', theColor*0.5);
    end


    % Render RGC RFs as semitransparent contours
    cMap = brewermap(1024, '*greys');
    rgcsNum = numel(RGCRFinputs);

    for iRGC = 1:rgcsNum 
        % Retrieve indices of input cones and their weights
        theConeInputs = RGCRFinputs{iRGC};
        if (isempty(theConeInputs))
            %fprintf(2,'Skip rendering RGC %d which has 0 inputs.', iRGC);
            continue;
        end

        if (isempty(RGCRFweights))
            theConeWeights = ones(1,numel(theConeInputs));
        else
            theConeWeights = RGCRFweights{iRGC};
        end

        %fprintf('Rendering RGC with %d cone inputs\n', numel(theConeWeights));

        % Generate RF outline
        xSupport = linspace(...
            min(allConePositions(theConeInputs,1))-min(allConeSpacings(theConeInputs)), ...
            max(allConePositions(theConeInputs,1))+min(allConeSpacings(theConeInputs)), ...
            64);
        ySupport = linspace(...
            min(allConePositions(theConeInputs,2))-min(allConeSpacings(theConeInputs)), ...
            max(allConePositions(theConeInputs,2))+min(allConeSpacings(theConeInputs)), ...
            64);

        [theRGCRFSpatialSupportXY, theRGCRF, theInputConeLineSpreadFunctionsXY] = ...
            RGCRFconnector.rfFromConeInputs(...
                allConePositions(theConeInputs,:), ...
                allConeSpacings(theConeInputs), ...
                theConeWeights, ...
                'xSupport', xSupport, ...
                'ySupport', ySupport);

        % Plot contour of RF outline
        zLevels(1) = 0.05*min(theConeWeights);
        zLevels(2) = max(theConeWeights);
        faceAlpha = 0.2;
        lineStyle = '-';
        lineWidth = 0.5;
        transparentContourPlot(ax, theRGCRFSpatialSupportXY, theRGCRF, ...
            zLevels, faceAlpha, cMap, lineStyle, lineWidth);

        if (superimposeConeInputWiring)
            if (numel(theConeInputs)>1)
                centroid = RGCRFconnector.centroidsFromConeInputs({theConeInputs}, {theConeWeights}, allConePositions);
                if (displayRGCID)
                    text(ax, centroid(1), centroid(2), sprintf('%d', iRGC), 'Color', [0 1 0], 'FontSize', 12, 'BackgroundColor', [0 0 0]);
                end
                for iInputCone = 1:numel(theConeInputs)
                    xx = [centroid(1) allConePositions(theConeInputs(iInputCone),1)];
                    yy = [centroid(2) allConePositions(theConeInputs(iInputCone),2)];
                    plot(ax, xx, yy, 'k-', 'LineWidth', theConeWeights(iInputCone)/max(theConeWeights)*3);
                end
            end
        end

        if (pauseAfterEachRGCisRendered)
            drawnow;
            pause
        end
        
    end

    % Finalize figure
    axis 'equal';

    dd = max(allConeSpacings);
    if (isempty(XLims))
        XLims = [min(allConePositions(:,1))-dd max(allConePositions(:,1))+dd];
    end
    if (isempty(YLims))
        YLims = [min(allConePositions(:,2))-dd max(allConePositions(:,2))+dd];
    end
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 16);
    title(ax, titleString);

    drawnow;
    if (~isempty(videoOBJ))
        videoOBJ.writeVideo(getframe(hFig));
    end

end

function [f,v] = facesAndVertices(positions, spacings, diskOutline)
    thetaSamples = size(diskOutline,1);
    rfsNum = size(positions, 1);
    X = zeros(rfsNum*thetaSamples,1);
    Y = X;
    f = zeros(rfsNum,thetaSamples);
    for iRF = 1:rfsNum
        ii = (iRF-1)*thetaSamples + (1:thetaSamples);
        f(iRF,:) = ii;
        X(ii,:) = positions(iRF,1) + spacings(iRF) * diskOutline(:,1);
        Y(ii,:) = positions(iRF,2) + spacings(iRF) * diskOutline(:,2);
    end
    v = [X(:) Y(:)];
end

function transparentContourPlot(axesHandle, spatialSupportXY, zData, ...
                                zLevels, faceAlpha, cmap, lineStyle, lineWidth)

    xSupport = squeeze(spatialSupportXY(:,1));
    ySupport = squeeze(spatialSupportXY(:,2));
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;
    hold(axesHandle, 'on');
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theLevelVerticesNum = C(2,startPoint);
        x = C(1,startPoint+(1:theLevelVerticesNum));
        y = C(2,startPoint+(1:theLevelVerticesNum));
        v = [x(:) y(:)];
        f = 1:numel(x);
        patch(axesHandle, 'Faces', f, 'Vertices', v, 'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(max([1 round(theLevel*size(cmap,1))]),:), ...
            'FaceAlpha', faceAlpha, ... 
            'EdgeAlpha', 1, ...
            'LineStyle', lineStyle, 'LineWidth', lineWidth);
        startPoint = startPoint + theLevelVerticesNum+1;
    end

end
