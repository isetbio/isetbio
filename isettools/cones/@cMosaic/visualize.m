function visualizationParams = visualize(obj, varargin)
% Visualize different aspects of a @cMosaic or its activation
%
% Syntax:
%   cm = cMosaic(); cm.visualize();
%
%   % Return the various settable params
%   pStruct = cm.visualize('params')
%
%   % Display the various settable params and info about them
%   cm.visualize('help');
%
%  Examples
%    cm.visualize();        - Visualize the cone mosaic
%    pStruct = cm.visualize('params') - Retrieve the visualization parameter structure
%
%  See also
%   Tutorials in tutorials/t_cones/cMosaic
%
    
    if ~isempty(varargin) && (isequal(varargin{1},'params') || isequal(varargin{1},'help'))
        visualizationParams = returnVisualizationParams(varargin{1});
        return;
    else
        visualizationParams = '';
    end

    % Parse input
    p = inputParser;
    p.addParameter('visualizationView', 'REVF', @(x)(ischar(x) && (ismember(x, {'REVF', 'retinal view'}))));
    p.addParameter('domain', 'degrees', @(x)(ischar(x) && (ismember(x, {'degrees', 'microns'}))));
    p.addParameter('domainVisualizationLimits', [], @(x)((isempty(x))||(numel(x)==4)));
    p.addParameter('domainVisualizationTicks', [], @(x)(isempty(x)||(isstruct(x)&&((isfield(x, 'x'))&&(isfield(x,'y'))))));
    
    p.addParameter('visualizedConeAperture', 'geometricArea', @(x)ismember(x, ...
        {'lightCollectingArea', 'geometricArea', 'coneSpacing', ...
        'lightCollectingAreaCharacteristicDiameter', 'lightCollectingArea2sigma', 'lightCollectingArea4sigma', 'lightCollectingArea5sigma', 'lightCollectingArea6sigma'}));
    p.addParameter('visualizedConeApertureThetaSamples', [], @(x)(isempty(x) || isscalar(x)));
    
    p.addParameter('visualizeCones', true, @islogical);
    p.addParameter('labelCones', true, @islogical);
    p.addParameter('labelConesInActivationMap', false, @islogical);
    p.addParameter('conesAlpha', 1.0, @isscalar);
    p.addParameter('conesEdgeAlpha', 1.0, @isscalar);
    p.addParameter('labelConesWithIndices', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('outlinedConesWithIndices', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('densityContourOverlay', false, @islogical);
    p.addParameter('densityContourLevels', [], @isnumeric);
    p.addParameter('densityContourLevelLabelsDisplay', false, @islogical);
    p.addParameter('densityColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    
    p.addParameter('withSuperimposedOpticalImage', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('superimposedOIAlpha', 0.7, @isnumeric);
    p.addParameter('withSuperimposedPSF', [], @(x)(isempty(x) || isstruct(x)));
    
    p.addParameter('activation', []);
    p.addParameter('activationRange', [],@(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('activationColorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.addParameter('verticalDensityColorBar', false, @islogical);
    p.addParameter('horizontalActivationColorBar', false, @islogical);
    p.addParameter('verticalActivationColorBar', false, @islogical);
    p.addParameter('horizontalActivationColorBarInside', false, @islogical);
    p.addParameter('verticalActivationColorBarInside', false, @islogical);
    p.addParameter('colorBarTickLabelPostFix', '', @ischar);
    p.addParameter('colorbarTickLabelColor',  [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    
    p.addParameter('horizontalActivationSliceEccentricity', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('verticalActivationSliceEccentricity', [], @(x)((isempty(x))||(isscalar(x))));
    
    p.addParameter('crossHairsOnMosaicCenter', false, @islogical);
    p.addParameter('crossHairsAtPosition', [], @(x)((isempty(x))||(numel(x)==2)));
    p.addParameter('crossHairsOnFovea', false, @islogical);
    p.addParameter('crossHairsOnOpticalImageCenter', false, @islogical);
    p.addParameter('crossHairsColor', [], @(x)(isempty(x)||((isvector(x))&&(numel(x) == 3))));
    
    p.addParameter('displayedEyeMovementData', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('currentEMposition', [], @(x)(isempty(x)||(numel(x)==2)));
    
    p.addParameter('labelRetinalMeridians', false, @islogical);
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('clearAxesBeforeDrawing', true, @islogical);
    p.addParameter('fontSize', 16, @isscalar);
    p.addParameter('fontAngle', 'normal', @(x)(ismember(lower(x), {'normal', 'italic'})));
    p.addParameter('colorbarFontSize', 16, @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('backgroundColor', [], @(x)( (ischar(x)&&((strcmp(x,'none'))||(strcmp(x,'mean of color map'))) ) || isempty(x) || ((isvector(x))&&(numel(x) == 3))));
    p.addParameter('plotTitle', '', @(x)(isempty(x) || ischar(x) || islogical(x)));
    p.addParameter('plotTitleColor', [0 0 0], @isnumeric);
    p.addParameter('plotTitleFontSize', 16, @isscalar);
    p.addParameter('textDisplay', '',@(x)(isempty(x) || ischar(x)));
    p.addParameter('textDisplayColor', [], @isnumeric);
    p.parse(varargin{:});
    
    % visualizationView = p.Results.visualizationView;
    domain = p.Results.domain;
    domainVisualizationLimits = p.Results.domainVisualizationLimits;
    domainVisualizationTicks = p.Results.domainVisualizationTicks;
    visualizedConeAperture = p.Results.visualizedConeAperture;
    visualizedConeApertureThetaSamples = p.Results.visualizedConeApertureThetaSamples;
    figureHandle = p.Results.figureHandle;
    axesHandle = p.Results.axesHandle;
    verticalDensityColorBar = p.Results.verticalDensityColorBar;
    densityContourOverlay = p.Results.densityContourOverlay;
    densityContourLevels = p.Results.densityContourLevels;
    densityContourLevelLabelsDisplay = p.Results.densityContourLevelLabelsDisplay;
    densityColorMap = p.Results.densityColorMap;
    superimposedOpticalImage = p.Results.withSuperimposedOpticalImage;
    superimposedPSF = p.Results.withSuperimposedPSF;
    activation = p.Results.activation;
    activationRange = p.Results.activationRange;
    currentEMposition = p.Results.currentEMposition;
    crossHairsOnMosaicCenter = p.Results.crossHairsOnMosaicCenter;
    crossHairsOnOpticalImageCenter = p.Results.crossHairsOnOpticalImageCenter;
    crossHairsAtPosition = p.Results.crossHairsAtPosition;
    visualizeCones = p.Results.visualizeCones;
    labelCones = p.Results.labelCones;
    labelConesInActivationMap = p.Results.labelConesInActivationMap;
    faceAlphaCones = p.Results.conesAlpha;
    edgeAlphaCones = p.Results.conesEdgeAlpha;
    labelConesWithIndices = p.Results.labelConesWithIndices;
    outlinedConesWithIndices = p.Results.outlinedConesWithIndices;
    labelRetinalMeridians = p.Results.labelRetinalMeridians;
    crossHairsOnFovea = p.Results.crossHairsOnFovea;
    crossHairsColor = p.Results.crossHairsColor;
    noXlabel = p.Results.noXLabel;
    noYlabel = p.Results.noYLabel;
    displayedEyeMovementData = p.Results.displayedEyeMovementData;
    fontSize = p.Results.fontSize;
    fontAngle = p.Results.fontAngle;
    plotTitleFontSize = p.Results.plotTitleFontSize;
    colorbarFontSize = p.Results.colorbarFontSize;
    cMap = p.Results.activationColorMap;
    verticalColorBar = p.Results.verticalActivationColorBar;
    colorbarTickLabelColor = p.Results.colorbarTickLabelColor;
    horizontalColorBar = p.Results.horizontalActivationColorBar;
    verticalColorBarInside = p.Results.verticalActivationColorBarInside;
    horizontalColorBarInside = p.Results.horizontalActivationColorBarInside;
    colorBarTickLabelPostFix = p.Results.colorBarTickLabelPostFix;
    horizontalActivationSliceEccentricity = p.Results.horizontalActivationSliceEccentricity;
    verticalActivationSliceEccentricity = p.Results.verticalActivationSliceEccentricity;
    backgroundColor = p.Results.backgroundColor;
    plotTitle = p.Results.plotTitle;
    plotTitleColor = p.Results.plotTitleColor;
    textDisplay = p.Results.textDisplay;
    textDisplayColor = p.Results.textDisplayColor;
    clearAxesBeforeDrawing = p.Results.clearAxesBeforeDrawing;

    if (~isempty(activation))
        labelCones = false;
        if (labelConesInActivationMap)
            labelCones = true;
        end
    end

    if (isempty(backgroundColor))
        backgroundColor = [0.7 0.7 0.7];
    end
    
    if (isempty(colorbarTickLabelColor))
        colorbarTickLabelColor = [1 0.5 0];
    end
    
    if (~isempty(labelConesWithIndices))
        labelCones = false;
    end
    
    % Determine what eye movement data have to be displayed
    if (isstruct(displayedEyeMovementData))
        if (ischar(displayedEyeMovementData.trial))
            switch(displayedEyeMovementData.trial)
                case 'all'
                    displayedTrials = 1:size(obj.fixEMobj.emPosArcMin,1);
                otherwise
                    error('unknown ''displayedEyeMovementData.trial'': ''%s''.', displayedEyeMovementData.trial);
            end
        else
            displayedTrials = displayedEyeMovementData.trial;
        end
        if (ischar(displayedEyeMovementData.timePoints))
            switch(displayedEyeMovementData.timePoints)
                case 'all'
                    displayedTimePoints = 1:size(obj.fixEMobj.emPosArcMin,2);
                otherwise
                    error('unknown ''displayedEyeMovementData.samples'': ''%s''.', displayedEyeMovementData.samples);
            end
        else
            displayedTimePoints = displayedEyeMovementData.timePoints;
        end
    end
    
    % Determine displayed domain (degs or microns)
    switch (domain)
        case 'degrees'
            rfPositions = obj.coneRFpositionsDegs;
            rfSpacings = obj.coneRFspacingsDegs;
            rfApertureDiameters = obj.coneApertureDiametersDegs;
            rfDiameters = rfApertureDiameters/obj.coneApertureToDiameterRatio;
            
            rfProximityThreshold = 1/270;
            if (isstruct(displayedEyeMovementData))
                emPath = -1/60*obj.fixEMobj.emPosArcMin(displayedTrials,displayedTimePoints,:);
            else
                emPath = [];
            end
        case 'microns'
            rfPositions = obj.coneRFpositionsMicrons;
            rfSpacings = obj.coneRFspacingsMicrons;
            rfApertureDiameters = obj.coneApertureDiametersMicrons;
            rfDiameters = rfApertureDiameters/obj.coneApertureToDiameterRatio;
            
            rfProximityThreshold = 1;
            if (isstruct(displayedEyeMovementData))
                emPath = -obj.fixEMobj.emPosMicrons(displayedTrials,displayedTimePoints,:);
            else
                emPath = [];
            end
    end
    
    
    % Show the emPath on the center of the mosaic
    emPath = bsxfun(@plus, emPath, reshape(mean(rfPositions,1), [1 1 2]));
    
    % Determine X,Y limits
    if (isempty(domainVisualizationLimits))
        xRange(1) = min(rfPositions(:,1));
        xRange(2) = max(rfPositions(:,1));
        if (xRange(2) == xRange(1))
            switch (domain)
                case 'degrees'
                    xRange = xRange(1) + 0.02*[-1 1];
                case 'microns'
                    xRange = xRange(1) + 2*[-1 1];
            end
        end
        yRange(1) = min(rfPositions(:,2));
        yRange(2) = max(rfPositions(:,2));
        if (yRange(2) == yRange(1))
            switch (domain)
                case 'degrees'
                    yRange = yRange(1) + 0.02*[-1 1];
                case 'microns'
                    yRange = yRange(1) + 2*[-1 1];
            end
        end
        xx = xRange(2)-xRange(1);
        yy = yRange(2)-yRange(1);
        xRange(1) = xRange(1)-xx*0.02;
        xRange(2) = xRange(2)+xx*0.02;
        yRange(1) = yRange(1)-yy*0.02;
        yRange(2) = yRange(2)+yy*0.02;
    else
        xRange(1) = domainVisualizationLimits(1);
        xRange(2) = domainVisualizationLimits(2);
        yRange(1) = domainVisualizationLimits(3);
        yRange(2) = domainVisualizationLimits(4);
    end
    
    
    % Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
        axesHandle = subplot('Position', [0.09 0.07 0.85 0.90]);
    else
        if (isempty(axesHandle))
            figure(figureHandle);
            clf;
            set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
            axesHandle = subplot('Position', [0.09 0.07 0.85 0.90]);
        end

        if (clearAxesBeforeDrawing)
            cla(axesHandle);
        end
    end
    visualizationParams.figureHandle = figureHandle;  % Return these.
    visualizationParams.axesHandle   = axesHandle;
    
    % Number of cones
    conesNum = numel(rfSpacings);
    
    % Aperture shape (disk)
    if (isempty(visualizedConeApertureThetaSamples))
        if (conesNum < 100)
            deltaAngle = 5;
        elseif (conesNum < 500)
            deltaAngle = 10;
        elseif (conesNum < 1000)
            deltaAngle = 15;
        elseif (conesNum < 2000)
            deltaAngle = 20;
        elseif (conesNum < 40000)
            deltaAngle = 30;
        elseif (conesNum < 60000)
            deltaAngle = 45;
        else
            deltaAngle = 60;
        end
    else
        deltaAngle = 360/visualizedConeApertureThetaSamples;
    end
    
    % Generate cone aperture shape
    iTheta = (0:deltaAngle:360) / 180 * pi;
    coneApertureShape.x = cos(iTheta);
    coneApertureShape.y = sin(iTheta);
    
    if (clearAxesBeforeDrawing)
        cla(axesHandle);
    end
    hold(axesHandle, 'on');
    
    
    % Visualize cone apertures
    switch (visualizedConeAperture)
        case 'coneSpacing'
            visualizedApertures = rfSpacings;
            
        case 'geometricArea'
            visualizedApertures = rfDiameters;
            
        case 'lightCollectingArea'
            visualizedApertures = rfApertureDiameters;
            
        case 'lightCollectingAreaCharacteristicDiameter'
            if (isfield(obj.coneApertureModifiers, 'shape') && (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.coneApertureModifiers.sigma;
                visualizedApertures = 2*sqrt(2)*gaussianSigma*rfApertureDiameters;
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize characteristic radius. Visualizing the diameter\n');
                visualizedApertures = rfApertureDiameters;
            end
            
        case 'lightCollectingArea2sigma'
            if (isfield(obj.coneApertureModifiers, 'shape') && (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.coneApertureModifiers.sigma;
                visualizedApertures = 2*gaussianSigma*rfApertureDiameters;
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize 2xsigma. Visualizing the diameter\n');
                visualizedApertures = rfApertureDiameters;
            end
            
        case 'lightCollectingArea4sigma'
            if (isfield(obj.coneApertureModifiers, 'shape') && (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.coneApertureModifiers.sigma;
                visualizedApertures = 4*gaussianSigma*rfApertureDiameters;
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize 4xsigma. Visualizing the diameter\n');
                visualizedApertures = rfApertureDiameters;
            end
            
        case 'lightCollectingArea5sigma'
            if (isfield(obj.coneApertureModifiers, 'shape') && (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.coneApertureModifiers.sigma;
                visualizedApertures = 5*gaussianSigma*rfApertureDiameters;
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize 5xsigma. Visualizing the diameter\n');
                visualizedApertures = rfApertureDiameters;
            end
            
        case 'lightCollectingArea6sigma'
            if (isfield(obj.coneApertureModifiers, 'shape') && (strcmp(obj.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.coneApertureModifiers.sigma;
                visualizedApertures = 6*gaussianSigma*rfApertureDiameters;
            else
                visualizedApertures = rfApertureDiameters;
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize 6xsigma. Visualizing the diameter\n');
            end
            
        otherwise
            error('Unknown visualizedConeAperture (%s)', visualizedConeAperture);
    end
    
    
    if (~isempty(activation))
        if (isempty(activationRange))
            activationRange(1) = min(activation(:));
            activationRange(2) = max(activation(:));
        end
        if (activationRange(1) == activationRange(2))
            activationRange(1) = activationRange(1)-0.1;
            activationRange(2) = activationRange(2)+0.1;
        end
        
        
        activation = (activation - activationRange(1))/(activationRange(2)-activationRange(1));
        activation(activation<0) = 0;
        activation(activation>1) = 1;
        
        % Visualize activations
        faceAlpha = 1.0;
        edgeAlpha = 1.0;
        % Plot L-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.lConeIndices)*0.5, ...
            rfPositions(obj.lConeIndices,:), activation(obj.lConeIndices), [0 0 0], 0.1, faceAlpha, edgeAlpha);
        % Plot M-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.mConeIndices)*0.5, ...
            rfPositions(obj.mConeIndices,:), activation(obj.mConeIndices), [0 0 0], 0.1, faceAlpha, edgeAlpha);
        % Plot S-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.sConeIndices)*0.5, ...
            rfPositions(obj.sConeIndices,:), activation(obj.sConeIndices), [0 0 0], 0.1, faceAlpha, edgeAlpha);
        % Plot K-cone activations
        renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.kConeIndices)*0.5, ...
            rfPositions(obj.kConeIndices,:), activation(obj.kConeIndices), [0 0 0], 0.1, faceAlpha, edgeAlpha);
        
        if (~isempty(verticalActivationSliceEccentricity))
            d = abs(rfPositions(:,2)-verticalActivationSliceEccentricity);
            idx = find(d < rfProximityThreshold);
            activationSlice = squeeze(activation(idx));
            h = stem(rfPositions(idx,1), rfPositions(idx,2) + activationSlice*0.2*(yRange(2)-yRange(1)), ...
                'filled',  'g-', 'LineWidth', 1.5);
            h.BaseValue = verticalActivationSliceEccentricity;
        end
    end

        
    if (visualizeCones)
        % Visualize cone types
        if (densityContourOverlay)
            faceAlpha = 0.2;
        else
            if (strcmp(visualizedConeAperture, 'lightCollectingArea6sigma'))
                faceAlpha = 0.6;
            else
                faceAlpha = 0.9;
            end
        end

        lineWidth = 1.0;
        % Plot L-cones
        if (labelCones) || (~isempty(labelConesWithIndices))
            edgeColor = [0.1 0.1 0.1];
        else
            edgeColor = [0.5 0.5 0.5];
        end
        
        if (~isempty(activation))
            faceAlphaCones = 0.0;
            edgeColor = [1 0 0];
            lineWidth = 0.5;
        end

        if (labelCones)
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.lConeIndices)*0.5, ...
                rfPositions(obj.lConeIndices,:), 1/4*0.9, edgeColor, lineWidth, faceAlphaCones, edgeAlphaCones);
        elseif (~isempty(labelConesWithIndices))
            includedLconeIndices = intersect(obj.lConeIndices, labelConesWithIndices);
            excludedLconeIndices = setdiff(obj.lConeIndices, includedLconeIndices);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(includedLconeIndices)*0.5, ...
                rfPositions(includedLconeIndices,:), 1/5*0.9, [1 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(excludedLconeIndices)*0.5, ...
                rfPositions(excludedLconeIndices,:), 5/4*0.9, [0 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
        end
        
        % Plot M-cones
        if (labelCones) || (~isempty(labelConesWithIndices))
            edgeColor = [0.1 0.1 0.1];
        else
            edgeColor = [0.5 0.5 0.5];
        end
        
        if (~isempty(activation))
            faceAlphaCones = 0.0;
            edgeColor = [0 1 0];
            lineWidth = 0.5;
        end

        if (labelCones)
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.mConeIndices)*0.5, ...
                rfPositions(obj.mConeIndices,:), 2/4*0.9, edgeColor, lineWidth, faceAlphaCones, edgeAlphaCones);
        elseif (~isempty(labelConesWithIndices))
            includedMconeIndices = intersect(obj.mConeIndices, labelConesWithIndices);
            excludedMconeIndices = setdiff(obj.mConeIndices, includedMconeIndices);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(includedMconeIndices)*0.5, ...
                rfPositions(includedMconeIndices,:), 2/5*0.9, [0 1 0], lineWidth, faceAlphaCones, edgeAlphaCones);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(excludedMconeIndices)*0.5, ...
                rfPositions(excludedMconeIndices,:), 5/4*0.9, [0 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
        end
        
        
        
        % Plot S-cones
        if (labelCones)  || (~isempty(labelConesWithIndices))
            edgeColor = [0.1 0.1 0.1];
        else
            edgeColor = [0.5 0.5 0.5];
        end

        if (~isempty(activation))
            faceAlphaCones = 0.0;
            edgeColor = [0 0.5 1];
            lineWidth = 0.5;
        end

        if (labelCones)
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.sConeIndices)*0.5, ...
                rfPositions(obj.sConeIndices,:), 3/4*0.9, edgeColor, lineWidth, faceAlphaCones, edgeAlphaCones);
        elseif (~isempty(labelConesWithIndices))
            includedSconeIndices = intersect(obj.sConeIndices, labelConesWithIndices);
            excludedSconeIndices = setdiff(obj.sConeIndices, includedSconeIndices);

            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(includedSconeIndices)*0.5, ...
                rfPositions(includedSconeIndices,:), 3/5*0.9, [0 0 1], lineWidth, faceAlphaCones, edgeAlphaCones);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(excludedSconeIndices)*0.5, ...
                rfPositions(excludedSconeIndices,:), 5/4*0.9, [0 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
        end
        
        
        % Plot K-cones
        if (labelCones)  || (~isempty(labelConesWithIndices))
            edgeColor = [1 1 0];
        else
            edgeColor = [0.5 0.5 0.5];
        end
        if (labelCones)
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(obj.kConeIndices)*0.5, ...
                rfPositions(obj.kConeIndices,:), 4/4*0.9, [0 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
        elseif (~isempty(labelConesWithIndices))
            includedKconeIndices = intersect(obj.kConeIndices, labelConesWithIndices);
            excludedKconeIndices = setdiff(obj.kConeIndices, includedKconeIndices);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(includedKconeIndices)*0.5, ...
                rfPositions(includedKconeIndices,:), 4/5*0.9, edgeColor, lineWidth, faceAlphaCones, edgeAlphaCones);
            renderPatchArray(axesHandle, coneApertureShape, visualizedApertures(excludedKconeIndices)*0.5, ...
                rfPositions(excludedKconeIndices,:), 5/4*0.9, [0 0 0], lineWidth, faceAlphaCones, edgeAlphaCones);
        end

        if (~isempty(outlinedConesWithIndices))
            for iOutlinedCone = 1:numel(outlinedConesWithIndices)
                theOutlinedConePos = rfPositions(outlinedConesWithIndices(iOutlinedCone),:);
                theOutlinedConeRectWidth = rfDiameters(outlinedConesWithIndices(iOutlinedCone));
                xx = theOutlinedConePos(1) + 0.5*theOutlinedConeRectWidth*[-1  1 1 -1 -1];
                yy = theOutlinedConePos(2) + 0.5*theOutlinedConeRectWidth*[-1 -1 1  1 -1];
                plot(axesHandle, xx, yy, 'k-', 'LineWidth', 2.0);
                plot(axesHandle, xx, yy, 'y--', 'LineWidth', 1.0);
            end
        end


    end % visualizeCones
    
    if (densityContourOverlay)
        % Compute dense 2D map
        sampledPositions{1} = linspace(xRange(1), xRange(2), 24);
        sampledPositions{2} = linspace(yRange(1), yRange(2), 24);
        
        % Convert spacing to density
        if (strcmp(domain, 'microns'))
            % Convert to mm, so we report density in cones / mm^2
            density2DMap = cMosaic.densityMap(rfPositions, rfSpacings/1e3, sampledPositions);
        else
            density2DMap = cMosaic.densityMap(rfPositions, rfSpacings, sampledPositions);
        end
        
        [densityContourX,densityContourY] = meshgrid(sampledPositions{1}, sampledPositions{2});
        
        % Render contour map
        if (isempty(densityContourLevels))
            densityContourLevels = round(prctile(density2DMap(:), [1 5 15 30 50 70 85 95 99])/100)*100;
        end
        contourLabelSpacing = 4000;
        if (isempty(densityColorMap))
            [cH, hH] = contour(axesHandle, densityContourX, densityContourY, ...
                density2DMap, densityContourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
                'ShowText', densityContourLevelLabelsDisplay, 'LabelSpacing', contourLabelSpacing);
            clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [0 0 0], 'BackgroundColor', 'none');
        else
            [cH, hH] = contourf(axesHandle, densityContourX, densityContourY, ...
                density2DMap, densityContourLevels, 'LineColor', 'k', 'LineWidth', 2.0, ...
                'ShowText', densityContourLevelLabelsDisplay, 'LabelSpacing', contourLabelSpacing);
            clabel(cH,hH,'FontWeight','bold', 'FontSize', 16, ...
                'Color', [1 1 1], 'BackgroundColor', 'none');
        end
        
    end
        
    
    
    % Add crosshairs
    if (~isempty(crossHairsAtPosition))
        crossHairsColor = [0 0 0];
        % Crosshairs centered on the middle of the mosaic
        xx1 = [xRange(1) xRange(2)];
        yy1 = crossHairsAtPosition(2)*[1 1];
        xx2 = crossHairsAtPosition(1)*[1 1];
        yy2 = [yRange(1) yRange(2)];
        plot(axesHandle, xx1, yy1, '-', 'Color', crossHairsColor, 'LineWidth', 1.5);
        plot(axesHandle, xx2, yy2,  '-', 'Color', crossHairsColor,'LineWidth', 1.5); 
    end
    
    if (crossHairsOnMosaicCenter) || (crossHairsOnOpticalImageCenter) || (crossHairsOnFovea)
        
        if (isempty(crossHairsColor))
            if (isempty(activation))
                if (strcmp(backgroundColor, 'none'))
                    crossHairsColor = [0 0 0];
                else
                    crossHairsColor = 1-backgroundColor;
                end
            else
                crossHairsColor = [1 0 0];
            end
        end
        
        if (crossHairsOnMosaicCenter)
            % Crosshairs centered on the middle of the mosaic
            xx1 = [xRange(1) xRange(2)];
            yy1 = mean(yRange)*[1 1];
            xx2 = mean(xRange)*[1 1];
            yy2 = [yRange(1) yRange(2)];
        elseif (crossHairsOnFovea)
            % Crosshairs centered on [0 0]
            switch (domain)
                case 'degrees'
                    xx1 = 30*[-1 1];
                    yy2 = 30*[-1 1];
                case 'microns'
                    xx1 = 30*[-1 1]*300;
                    yy2 = 30*[-1 1]*300;
            end
            yy1 = [0 0];
            xx2 = [0 0];
        else
            % Crosshairs centered on the middle of the optical image, i.e.,
            % (0,0) or current em position, if that is passed in
            xx1 = [xRange(1) xRange(2)];
            if (~isempty(currentEMposition))
                yy1 = -[1 1] * currentEMposition(2);
                xx2 = -[1 1] * currentEMposition(1);
            else
                yy1 = [0 0];
                xx2 = [0 0];
            end
            yy2 = [yRange(1) yRange(2)];
        end
        plot(axesHandle, xx1, yy1, '-', 'Color', crossHairsColor, 'LineWidth', 1.5);
        plot(axesHandle, xx2, yy2,  '-', 'Color', crossHairsColor,'LineWidth', 1.5);
    end
    
    % Superimpose eye movement path(s)
    if (~isempty(emPath))
        assert(size(emPath,3) == 2, sprintf('The third dimension of an emPath must be 2.'));
        nTrials = size(emPath,1);
        lineColors = brewermap(nTrials, 'Spectral');
        for emTrialIndex = 1:nTrials
            plot(axesHandle, emPath(emTrialIndex,:,1), emPath(emTrialIndex,:,2), ...
                'k-', 'LineWidth', 4.0);
        end
        for emTrialIndex = 1:nTrials
            plot(axesHandle, emPath(emTrialIndex,:,1), emPath(emTrialIndex,:,2), ...
                '-', 'LineWidth', 2, 'Color', squeeze(lineColors(emTrialIndex,:)));
        end
    end
    
    % Superimpose optical image
    if (~isempty(superimposedOpticalImage))
        superimposeTheOpticalImage(obj, axesHandle, domain, superimposedOpticalImage,p.Results.superimposedOIAlpha);
    end

    % Superimpose PSF
    if (~isempty(superimposedPSF))
        superimposeThePSF(obj, axesHandle, domain, superimposedPSF);
    end

    hold(axesHandle, 'off');
    
    % Set appropriate colormap
    if (isempty(activation))
        % Colormap for visualization of cone types
        if (labelCones)
            cMap = [obj.lConeColor; obj.mConeColor; obj.sConeColor; obj.kConeColor];
        elseif (~isempty(labelConesWithIndices))
            cMap = [obj.lConeColor; obj.mConeColor; obj.sConeColor; obj.kConeColor; [0.5 0.5 0.5]];
        else
            cMap = 0.4*ones(4,3);
        end
    else
        % Colormap for visualization of activity
        if (isempty(cMap))
            cMap = gray(numel(obj.coneRFspacingsDegs)); %brewermap(numel(obj.coneRFspacingsDegs), '*greys');
        end
        
        if (ischar(backgroundColor) && strcmp(backgroundColor, 'mean of color map'))
            midRow = round(size(cMap,1)/2);
            backgroundColor = squeeze(cMap(midRow,:));
        elseif (isempty(backgroundColor))
            backgroundColor = squeeze(cMap(1,:));
        end
    end
    
    
    
    if (~isempty(activation))
        if (verticalColorBar) || (horizontalColorBar) || (verticalColorBarInside) || (horizontalColorBarInside)
            colorBarTicks = [0.00 0.25 0.5 0.75 1.0];
            colorBarTickLabels = cell(1, numel(colorBarTicks));
            colorBarTickLevels = activationRange(1) + (activationRange(2)-activationRange(1)) * colorBarTicks;
            
            for k = 1:numel(colorBarTicks)
                if (max(abs(colorBarTickLevels)) >= 10)
                    colorBarTickLabels{k} = sprintf('%2.0f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 1)
                    colorBarTickLabels{k} = sprintf('%2.1f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                elseif (max(abs(colorBarTickLevels)) >= 0.1)
                    colorBarTickLabels{k} = sprintf('%2.2f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                else
                    colorBarTickLabels{k} = sprintf('%2.3f %s', colorBarTickLevels(k), colorBarTickLabelPostFix);
                end
            end
            
            if (isempty(colorbarFontSize))
                colorbarFontSize = fontSize/2;
            end
    
            if (verticalColorBar)
                colorbar(axesHandle, 'eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (verticalColorBarInside)
                colorbar(axesHandle, 'east', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            elseif (horizontalColorBar)
                colorbar(axesHandle,'northOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor);
            elseif (horizontalColorBarInside)
                colorbar(axesHandle,'north', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels, ...
                    'Color', colorbarTickLabelColor,  'FontWeight', 'Bold', 'FontSize', colorbarFontSize, 'FontName', 'Spot mono');
            end
        else
            colorbar(axesHandle, 'off');
        end
    else
        if (verticalDensityColorBar)
            colorBarTicks = densityContourLevels;
            colorBarTickLabels = colorBarTicks;
            colorbar(axesHandle, 'eastOutside', 'Ticks', colorBarTicks, 'TickLabels', colorBarTickLabels);
        else
            colorbar(axesHandle, 'off');
        end
    end
    
    
    
    % Finalize plot
    xtickangle(axesHandle, 0);
    set(axesHandle, 'Color', backgroundColor);
    axis(axesHandle, 'xy');
    axis(axesHandle, 'equal');
    set(axesHandle, 'XLim', xRange, 'YLim', yRange,'FontSize', fontSize, 'FontAngle', fontAngle);
    
    if (~isempty(densityColorMap)) && (densityContourOverlay)
        colormap(axesHandle, densityColorMap);
        set(axesHandle, 'CLim', [min(densityContourLevels(:)) max(densityContourLevels(:))]);
    else
        colormap(axesHandle, cMap);
        set(axesHandle, 'CLim', [0 1]);
    end
    
    
    if (isempty(domainVisualizationTicks))
        xo = (xRange(1)+xRange(2))/2;
        xx = xRange(2)-xRange(1);
        yo = (yRange(1)+yRange(2))/2;
        yy = yRange(2)-yRange(1);
        ticksX = xo + xx*0.5*[-0.75 0 0.75];
        ticksY = yo + yy*0.5*[-0.75 0 0.75];
        
        if (xx > 10)
            domainVisualizationTicks.x = round(ticksX);
        elseif (xx > 5)
            domainVisualizationTicks.x = round(ticksX*10)/10;
        elseif (xx > 1)
            domainVisualizationTicks.x = round(ticksX*100)/100;
        else
            domainVisualizationTicks.x = round(ticksX*1000)/1000;
        end

        if (yy > 10)
            domainVisualizationTicks.y = round(ticksY);
        elseif (yy > 5)
            domainVisualizationTicks.y = round(ticksY*10)/10;
        elseif (yy > 1)
            domainVisualizationTicks.y = round(ticksY*100)/100;
        else
            domainVisualizationTicks.y = round(ticksY*1000)/1000;
        end
        
    end
    
    set(axesHandle, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
    
    box(axesHandle, 'on');
    set(figureHandle, 'Color', [1 1 1]);
    
    switch (domain)
        case 'degrees'
            if (~noXlabel)
                if (labelRetinalMeridians)
                    if (strcmp(obj.whichEye, 'right eye'))
                        leftMeridianName = 'temporal retina';
                        rightMeridianName = 'nasal retina';
                    else
                        leftMeridianName = 'nasal retina';
                        rightMeridianName = 'temporal retina';
                    end
                    xlabel(axesHandle, sprintf('\\color{red}%s    \\color{black} space (degrees)    \\color[rgb]{0 0.7 0} %s', ...
                        leftMeridianName, rightMeridianName));
                else
                    xlabel(axesHandle, 'space (degrees)');
                end
            end
            if (~noYlabel)
                if (labelRetinalMeridians)
                    ylabel(axesHandle, sprintf('%s  < = = = = = |     space (degrees)    | = = = = =  > %s', ...
                        'superior retina', 'inferior retina'));
                else
                    ylabel(axesHandle, 'space (degrees)');
                end
            end
            minTickIncrement = min([min(abs(diff(domainVisualizationTicks.x))) min(abs(diff(domainVisualizationTicks.y)))]);
            
            if (minTickIncrement >= 0.1)
                set(axesHandle, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
                    'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
            else
                set(axesHandle, 'XTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.x), ...
                    'YTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.y));
            end
            
        case 'microns'
            if (~noXlabel)
                if (labelRetinalMeridians)
                    if (strcmp(obj.whichEye, 'right eye'))
                        leftMeridianName = '(temporal)';
                        rightMeridianName = '(nasal)';
                    else
                        leftMeridianName = '(nasal)';
                        rightMeridianName = '(temporal)';
                    end
                    xlabel(axesHandle, sprintf('\\color{red}%s    \\color{black} space (microns)    \\color[rgb]{0 0.7 0} %s', ...
                        leftMeridianName, rightMeridianName));
                else
                    xlabel(axesHandle, 'space (microns)');
                end
            end
            if (~noYlabel)
                if (labelRetinalMeridians)
                    upperMeridianName = '(inferior)';
                    lowerMeridianName = '(superior)';
                    ylabel(axesHandle, sprintf('\\color{blue}%s    \\color{black} space (microns)    \\color[rgb]{0.6 0.6 0.4} %s', ...
                        lowerMeridianName, upperMeridianName));
                else
                    ylabel(axesHandle, 'space (microns)');
                end
            end
            minTickIncrement = min([min(abs(diff(domainVisualizationTicks.x))) min(abs(diff(domainVisualizationTicks.y)))]);
            if (minTickIncrement >= 1)
                set(axesHandle, 'XTickLabel', sprintf('%1.0f\n', domainVisualizationTicks.x), ...
                    'YTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.y));
            else
                set(axesHandle, 'XTickLabel', sprintf('%1.1f\n', domainVisualizationTicks.x), ...
                    'YTickLabel', sprintf('%1.2f\n', domainVisualizationTicks.y));
            end
            
    end
    
    % User can set plotTitle to false.  Then we skip.  Used in
    % cMosaic.plot
    if plotTitle    
        if (isempty(plotTitle))
            if (numel(obj.coneDensities) == 4)
                title(axesHandle,sprintf('L (%2.1f%%), M (%2.1f%%), S (%2.1f%%), K (%2.1f%%), N = %d', ...
                    100*obj.coneDensities(1), ...
                    100*obj.coneDensities(2), ...
                    100*obj.coneDensities(3), ...
                    100*obj.coneDensities(4), ...
                    conesNum), 'Color', plotTitleColor, 'FontSize', plotTitleFontSize);
            else
                title(axesHandle,sprintf('L (%2.1f%%), M (%2.1f%%), S (%2.1f%%), N = %d', ...
                    100*obj.coneDensities(1), ...
                    100*obj.coneDensities(2), ...
                    100*obj.coneDensities(3), ...
                    conesNum), 'Color', plotTitleColor, 'FontSize', plotTitleFontSize);
            end
        else
            title(axesHandle,plotTitle, 'Color', plotTitleColor, 'FontSize', plotTitleFontSize);
        end
    end
    
    if (~isempty(textDisplay))
        dx = 0.45*(xRange(2)-xRange(1));
        dy = 0.04*(yRange(2)-yRange(1));
        if (isempty(textDisplayColor))
            textDisplayColor = 1-backgroundColor;
        end
        text(axesHandle, xRange(1)+dx, yRange(1)+dy, textDisplay, ...
            'FontSize', 16, 'Color', textDisplayColor, 'BackgroundColor', backgroundColor);
    end
    
    drawnow;
end

function superimposeThePSF(obj, axesHandle, visualizationDomain, thePSFData)

    xSupport = thePSFData.supportXdegs + obj.eccentricityDegs(1);
    ySupport = thePSFData.supportYdegs + obj.eccentricityDegs(2);

    
    if (strcmp(visualizationDomain, 'microns'))
        if (isfield(thePSFData, 'supportXmicrons'))
            xSupport  = thePSFData.supportXmicrons + obj.eccentricityMicrons(1);
            ySupport  = thePSFData.supportYmicrons + obj.eccentricityMicrons(2);
        else
             % Convert spatial support in microns to degs
             xSupport  = obj.distanceDegreesToDistanceMicronsForCmosaic(xSupport);
             ySupport  = obj.distanceDegreesToDistanceMicronsForCmosaic(ySupport);
        end
    end

    cmap = brewermap(1024,'blues');
    alpha = 0.75;
    contourLineColor = [0.2 0.2 0.2];

    cMosaic.semiTransparentContourPlot(axesHandle, ...
            xSupport, ySupport, ...
            thePSFData.data/max(thePSFData.data(:)), ...
            0.05:0.15:0.95, cmap, alpha, contourLineColor, ...
            'lineWidth', 1.5);
end

function superimposeTheOpticalImage(obj, axesHandle, visualizationDomain, theOI, superimposedOIAlpha)

    % Obtain spatial support in microns
    spatialSupportMeters = oiGet(theOI, 'spatial support');
    xSupport = squeeze(spatialSupportMeters(1,1:end,1)) * 1e6;
    ySupport = squeeze(spatialSupportMeters(1:end,1,2)) * 1e6;

    if (strcmp(visualizationDomain, 'degrees'))
         % Convert spatial support in microns to degs
         xSupport  = obj.distanceMicronsToDistanceDegreesForCmosaic(xSupport);
         ySupport  = obj.distanceMicronsToDistanceDegreesForCmosaic(ySupport);
    end

    % Translate to position
    if ((ischar(obj.opticalImagePositionDegs))&&(strcmp(obj.opticalImagePositionDegs, 'mosaic-centered')))
        if (strcmp(visualizationDomain, 'degrees'))
            xSupport = xSupport + obj.eccentricityDegs(1);
            ySupport = ySupport + obj.eccentricityDegs(2);
        else
            xSupport = xSupport + obj.eccentricityMicrons(1);
            ySupport = ySupport + obj.eccentricityMicrons(2);
        end
    else
        % Translate to custom position
        x0 = obj.opticalImagePositionDegs(1);
        y0 = obj.opticalImagePositionDegs(2);
        if (strcmp(visualizationDomain, 'microns'))
            % Convert degs to microns
            x0  =  obj.distanceMicronsToDistanceDegreesForCmosaic(x0);
            x0  =  obj.distanceMicronsToDistanceDegreesForCmosaic(y0);
        end
        xSupport = xSupport + x0;
        ySupport = ySupport + y0;
    end

    % Get the rgb value
    theOpticalImageRGB = flipud(oiGet(theOI, 'rgb'));
    % Superimpose the optical image
    imPlot = image(axesHandle, xSupport, ySupport, theOpticalImageRGB);
    % Make it semi-transparent
    imPlot.AlphaData = superimposedOIAlpha;
end


function renderPatchArray(axesHandle, apertureShape, apertureRadii, rfCoords, ...
    faceColors, edgeColor, lineWidth, faceAlpha, edgeAlpha)

    conesNum = numel(apertureRadii);
    if (conesNum == 0)
        return;
    end

    verticesPerCone = numel(apertureShape.x);
    
    verticesList = zeros(verticesPerCone * conesNum, 2);
    facesList = [];
    
    if (numel(faceColors) == 1)
        colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
    else
        colors = [];
    end
    
    for coneIndex = 1:conesNum
        idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
        verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + rfCoords(coneIndex,1);
        verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + rfCoords(coneIndex,2);
        if ((numel(faceColors) == conesNum)&& (conesNum > 1))
            colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
        end
        facesList = cat(1, facesList, idx);
    end
    
    S.Vertices = verticesList;
    S.Faces = facesList;
    S.FaceVertexCData = colors;
    S.FaceColor = 'flat';
    S.EdgeColor = edgeColor;
    S.FaceAlpha = faceAlpha;
    S.EdgeAlpha = edgeAlpha;
    S.LineWidth = lineWidth;
    patch(S, 'Parent', axesHandle);
end


function params = returnVisualizationParams(mode)
    % User wants to return a struct with a list of parameters
    % This can be set and passed in as the varargin.
    
    visualizationParamsStruct.domain = struct(...
        'default', 'degrees', ...
        'docString', 'Choose between {''degrees'', ''microns''}');
    
    visualizationParamsStruct.domainVisualizationLimits = struct(...
        'default', [], ...
        'docString', 'Limits for visualization. Either [], or a 4-element vector, [xMin xMax yMin yMax]' ...
        );
    
    visualizationParamsStruct.domainVisualizationTicks = struct(...
        'default', [], ...
        'docString', 'Ticks for visualization. Either [], or a struct with the following content, struct(''x'', -10:1:10, ''y'', -5:1:5)' ...
        );
    
    visualizationParamsStruct.visualizedConeAperture = struct(...
        'default', 'geometricArea', ...
        'docStringA', 'For pillbox apertures,  choose between {''geometricArea'', ''coneSpacing'', ''lightCollectingArea''}', ...
        'docStringB', 'For Gaussian apertures, choose between {''geometricArea'', ''coneSpacing'', ''lightCollectingAreaCharacteristicDiameter'', ''lightCollectingArea2sigma'', ''lightCollectingArea4sigma'', ''lightCollectingArea5sigma'', ''lightCollectingArea6sigma''}' ...
        );
    
    visualizationParamsStruct.visualizedConeApertureThetaSamples = struct(...
        'default', [], ...
        'docString', 'Number of angular samples to represent the visualized cone aperture, e.g., 6 for a hexagonal-shaped aperture' ...
        );
    
    visualizationParamsStruct.visualizeCones = struct(...
        'default', true, ...
        'docString', 'Flag indicating whether to visual cones. You can set it to false when visualizing the cone density' ...
        );
    
    visualizationParamsStruct.labelCones = struct(...
        'default', true, ...
        'docString', 'Flag indicating whether to color-code the visualized cones according to the their type: L, M, S etc' ...
        );
    
    visualizationParamsStruct.labelConesWithIndices = struct(...
        'default', [], ...
        'docString', 'Either [], in which case all cones are labeled, or a list of indices specifying the cones to be labeled according to their type' ...
        );
    
    visualizationParamsStruct.densityContourOverlay  = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to superimpose a contour map of cone density on top of the cone mosaic' ...
        );
    
    visualizationParamsStruct.densityContourLevels = struct(...
        'default', [], ...
        'docString', 'Either [] or a vector of levels for which to draw iso-density contours' ...
        );
    
    visualizationParamsStruct.densityContourLevelLabelsDisplay = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to label the iso-density contours' ...
        );
    
    visualizationParamsStruct.densityColorMap = struct(...
        'default', [], ...
        'docString', 'Either [] or a [nx3] matrix of RGB values for the density map' ...
        );
    
    visualizationParamsStruct.activation = struct(...
        'default', [], ...
        'docString', 'Either [] or a [1xnCones] matrix of activation values returned by cMosaic.compute()' ...
        );
    
    visualizationParamsStruct.activationRange = struct(...
        'default', [], ...
        'docString', 'Either [] or a [1x2] vector of the visualized minimum and maximum activation values' ...
        );
    
    visualizationParamsStruct.activationColorMap  = struct(...
        'default', [], ...
        'docString', 'Either [] or a [nx3] matrix of RGB values for the activation map' ...
        );
    
    visualizationParamsStruct.verticalDensityColorBar = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add a vertically-oriented colorbar showing the color-coding of the density map' ...
        );
    
    visualizationParamsStruct.verticalActivationColorBar = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add a vertically-oriented colorbar showing the color-coding of the activation map' ...
        );
    visualizationParamsStruct.horizontalActivationColorBar = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add a horizontally-oriented colorbar showing the color-coding of the activation map' ...
        );
    
    visualizationParamsStruct.verticalActivationColorBarInside = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add a vertically-oriented colorbar showing the color-coding of the activation map. This colorbar is placed inside the density map.' ...
        );
    
    visualizationParamsStruct.horizontalActivationColorBarInside = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add a horizontally-oriented colorbar showing the color-coding of the activation map. This colorbar is placed inside the density map.' ...
        );
    
    visualizationParamsStruct.colorBarTickLabelPostFix = struct(...
        'default', '', ...
        'docString', 'Either an empty string, or a string to append to the labels of the colorbar' ...
        );
    
    visualizationParamsStruct.colorbarTickLabelColor = struct(...
        'default', [.9 .6 0.2], ...
        'docString', 'Either [], or a [1x3] vector of RGB values for the color of the colorbar ticks' ...
        );
    
    visualizationParamsStruct.horizontalActivationSliceEccentricity = struct(...
        'default', [], ...
        'docString', 'Either empty or a scalar specifying the y-eccentricity through which to depict the activation of cones along the x-axis' ...
        );
    
    visualizationParamsStruct.verticalActivationSliceEccentricity = struct(...
        'default', [], ...
        'docString', 'Either empty or a scalar specifying the x-eccentricity through which to depict the activation of cones along the y-axis' ...
        );
    
    visualizationParamsStruct.crossHairsOnMosaicCenter = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add cross-hairs located at the mosaic center' ...
        );
    
    visualizationParamsStruct.crossHairsOnFovea = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add cross-hairs located at the fovea' ...
        );
    
    visualizationParamsStruct.crossHairsOnOpticalImageCenter = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to add cross-hairs located at the center of the optical image' ...
        );
    
    visualizationParamsStruct.crossHairsColor = struct(...
        'default', [], ...
        'docString', 'Either [], or a [1x3] vector of RGB values for the color of the cross-hairs' ...
        );
    
    visualizationParamsStruct.displayedEyeMovementData = struct(...
        'default', [], ...
        'docStringA', 'Either [], or a struct specifying which eye movement data (trial index & time points) to superimpose.', ...
        'docStringB', 'The struct should have the following format: struct(''trial'', 1:nTrials, ''timePoints'', timePoint1:timePoint2).', ...
        'docStringC', 'See t_cMosaicSingleEyeMovementPath.m for usage'...
        );
    
    visualizationParamsStruct.currentEMposition = struct(...
        'default', [], ...
        'docStringA', 'Either [], or a [1x2] vector of the (x,y) coordinates of the current eye position', ...
        'docStringB', 'See t_cMosaicSingleEyeMovementPath.m for usage'...
        );
    
    visualizationParamsStruct.labelRetinalMeridians = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to label the corresponding retinal meridians.' ...
        );
    
    visualizationParamsStruct.noXLabel = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to avoid labeling the x-axis.' ...
        );
    
    visualizationParamsStruct.noYLabel = struct(...
        'default', false, ...
        'docString', 'Flag indicating whether to avoid labeling the y-axis.' ...
        );
    
    visualizationParamsStruct.figureHandle = struct(...
        'default', [], ...
        'docString', 'Either [], or the handle to an existing figure on which to render the mosaic' ...
        );
    
    visualizationParamsStruct.axesHandle = struct(...
        'default', [], ...
        'docString', 'Either [], or the handle to existing axes on which to render the mosaic' ...
        );
    
    visualizationParamsStruct.fontSize = struct(...
        'default', 16, ...
        'docString', 'The font size for the figure' ...
        );
    
    visualizationParamsStruct.backgroundColor = struct(...
        'default', [], ...
        'docString', 'Either [], ''none'', or a [1x3] vector of RGB values for the axes background color' ...
        );
    
    visualizationParamsStruct.plotTitle = struct(...
        'default', [], ...
        'docString', 'Either an empty string, or a string specifying the title displayed above the mosaic' ...
        );
    
    visualizationParamsStruct.textDisplay = struct(...
        'default', [], ...
        'docString', 'Either an empty string, or a string specifying the text to be displayed at the bottom of the mosaic' ...
        );
    
    visualizationParamsStruct.textDisplayColor = struct(...
        'default', [], ...
        'docString', 'Either [], or a [1x3] vector of RGB values for the displayed text' ...
        );
    
    
    switch (mode)
        case 'params'
            % Return the params struct
            fNames = fieldnames(visualizationParamsStruct);
            for k = 1:numel(fNames)
                params.(fNames{k}) = visualizationParamsStruct.(fNames{k}).default;
            end
            
        case 'help'
            % Display the params struct along with information for each param
            visualizeParamsStructTree(visualizationParamsStruct, 'visualizationParams');
            params = '';
    end
end
