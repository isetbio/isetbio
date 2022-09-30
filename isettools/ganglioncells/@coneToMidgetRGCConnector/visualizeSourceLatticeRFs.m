function visualizeSourceLatticeRFs(obj, ax, coneOutline, varargin)
    p = inputParser;
    p.addParameter('identifiedConeIndicesSetA', [], @(x)(isempty(x)||isnumeric(x)));
    p.addParameter('identifiedConeIndicesSetB', [], @(x)(isempty(x)||isnumeric(x)));
    p.parse(varargin{:});
    identifiedConeIndicesSetA = p.Results.identifiedConeIndicesSetA;
    identifiedConeIndicesSetB = p.Results.identifiedConeIndicesSetB;

    % We have metadata, so plot the types of cones
    if (isfield(obj.sourceLattice, 'metaData')) && ...
           (isfield(obj.sourceLattice.metaData, 'coneTypes')) && ...
           (isfield(obj.sourceLattice.metaData, 'coneTypeIDs')) && ...
           (isfield(obj.sourceLattice.metaData, 'coneColors'))

        for iConeType = 1:numel(obj.sourceLattice.metaData.coneTypeIDs)
            % Plot each type separately
            idx = find(obj.sourceLattice.metaData.coneTypes == obj.sourceLattice.metaData.coneTypeIDs(iConeType));
            theColor = squeeze(obj.sourceLattice.metaData.coneColors(iConeType,:));
            
            % The positions
            coneRFpositionsMicrons = obj.sourceLattice.RFpositionsMicrons(idx,:);

            % The apertures
            coneRFspacingsMicrons = obj.sourceLattice.RFspacingsMicrons(idx);

            % Generate patch
            [f,v] = obj.facesAndVertices(coneRFpositionsMicrons, coneRFspacingsMicrons, coneOutline);
    
            if ((~isempty(identifiedConeIndicesSetA)) || (~isempty(identifiedConeIndicesSetB)))
                patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
            else
                patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'LineWidth', 1.0);
            end
            hold(ax, 'on')
        end
    else
        idx = 1:numel(obj.sourceLattice.RFspacingsMicrons);
        theColor = [0.3 0.3 0.3];
        
        % The positions
        coneRFpositionsMicrons = obj.sourceLattice.RFpositionsMicrons(idx,:);

        % The apertures
        coneRFaperturesMicrons = obj.sourceLattice.RFspacingsMicrons(idx);

        % Generate patch
        [f,v] = obj.facesAndVertices(coneRFpositionsMicrons, coneRFaperturesMicrons, coneOutline);

        if ((~isempty(identifiedConeIndicesSetA)) || (~isempty(identifiedConeIndicesSetB)))
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'FaceAlpha', 0.3, 'LineWidth', 1.0);
        else
            patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor*0.99, 'EdgeColor', theColor*0.7, 'LineWidth', 1.0);
        end
        hold(ax, 'on')

    end
end


