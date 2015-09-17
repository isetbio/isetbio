
        % Adding one layer
        function obj = addLayer(obj, cellType, ang)
            % Add a layer of ganglion cells to the rgcP structure
            % The properties of the different types are stored in the 
            % function layerSetCellType
            if notDefined('cellType')
                cellType = 'default';
                disp('Using "default" layer type.');
            end
            if notDefined('ang'), ang = []; end
            
            % Creating the layer
            layer = rgcLayer(obj);
            
            % Setting it to the wanted type
            layerSetCellType(layer,cellType,ang);
            
            % putting it into the layers field
            obj.layers{obj.get('nLayers') + 1} = layer;
        end