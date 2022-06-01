function visualizeConePositions(obj, ax, coneOutline, varargin)
% Visualize the cones of the input cone mosaic using a custom shape
% cone outline

    p = inputParser;
    p.addParameter('identifiedConeIndicesSetA', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('identifiedConeIndicesSetB', [], @(x)(isempty(x)||isnumeric(x)));
    p.parse(varargin{:});
    identifiedConeIndicesSetA = p.Results.identifiedConeIndicesSetA;
    identifiedConeIndicesSetB = p.Results.identifiedConeIndicesSetB;

    % Plot the cones
    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    coneColors = [obj.inputConeMosaic.lConeColor; obj.inputConeMosaic.mConeColor; obj.inputConeMosaic.sConeColor];
    
    for iConeType = 1:numel(coneTypes)
        idx = find(obj.inputConeMosaic.coneTypes == coneTypes(iConeType));
        allConeRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(idx,:);
        allConeRFspacingsMicrons = obj.inputConeMosaic.coneRFspacingsMicrons(idx);
        [f,v] = RGCconnector.facesAndVertices(allConeRFpositionsMicrons, allConeRFspacingsMicrons, coneOutline);
        theColor = squeeze(coneColors(iConeType,:));
        if ((~isempty(identifiedConeIndicesSetA)) || (~isempty(identifiedConeIndicesSetB)))
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor, 'EdgeColor', theColor*0.5, 'FaceAlpha', 0.3);
        else
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor, 'EdgeColor', theColor*0.5);
        end
        hold(ax, 'on')
    end

    % Plot identified cones using heavy lines
    if (~isempty(identifiedConeIndicesSetA))
        identifiedConeRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(identifiedConeIndicesSetA,:);
        identifiedConeRFspacingsMicrons = obj.inputConeMosaic.coneRFspacingsMicrons(identifiedConeIndicesSetA);
        [f,v] = RGCconnector.facesAndVertices(identifiedConeRFpositionsMicrons, identifiedConeRFspacingsMicrons, coneOutline);
        theColor = squeeze(coneColors(iConeType,:));
        patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', 'none', 'EdgeColor', theColor*0.5, 'LineWidth', 2.0);
    end

    if (~isempty(identifiedConeIndicesSetB))
        thetas = linspace(0,360,5);
        squareOutline = 0.5*[cosd(thetas); sind(thetas)]';

        identifiedConeRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(identifiedConeIndicesSetB,:);
        identifiedConeRFspacingsMicrons = obj.inputConeMosaic.coneRFspacingsMicrons(identifiedConeIndicesSetB);
        [f,v] = RGCconnector.facesAndVertices(identifiedConeRFpositionsMicrons, identifiedConeRFspacingsMicrons, squareOutline);
        theColor = squeeze(coneColors(iConeType,:));
        patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', 'none', 'EdgeColor', theColor*0.5, 'LineWidth', 2.0);

    end


end
