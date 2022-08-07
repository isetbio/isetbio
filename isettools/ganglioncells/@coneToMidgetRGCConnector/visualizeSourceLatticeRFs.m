function visualizeSourceLatticeRFs(obj, ax, coneOutline, varargin)
    p = inputParser;
    p.addParameter('identifiedConeIndicesSetA', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('identifiedConeIndicesSetB', [], @(x)(isempty(x)||isnumeric(x)));
    p.parse(varargin{:});
    identifiedConeIndicesSetA = p.Results.identifiedConeIndicesSetA;
    identifiedConeIndicesSetB = p.Results.identifiedConeIndicesSetB;


    if (~isfield(obj.sourceLattice, 'metaData'))
        % No meta data, so we cannot plot the types of cones
        return;
    end

    % We have metadata, so plot the types of cones
    for iConeType = 1:numel(obj.sourceLattice.metaData.coneTypeIDs)

        % Plot each type separately
        idx = find(obj.sourceLattice.metaData.coneTypes == obj.sourceLattice.metaData.coneTypeIDs(iConeType));
        
        % The positions
        coneRFpositionsMicrons = obj.sourceLattice.RFpositionsMicrons(idx,:);

        % The apertures
        coneRFaperturesMicrons = obj.sourceLattice.RFspacingsMicrons(idx);

        [f,v] = obj.facesAndVertices(coneRFpositionsMicrons, coneRFaperturesMicrons, coneOutline);
        theColor = squeeze(obj.sourceLattice.metaData.coneColors(iConeType,:));

        if ((~isempty(identifiedConeIndicesSetA)) || (~isempty(identifiedConeIndicesSetB)))
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
        else
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'LineWidth', 1.0);
        end
        hold(ax, 'on')

    end

    dataOut.coneRFaperturesMicrons = coneRFaperturesMicrons;
end
