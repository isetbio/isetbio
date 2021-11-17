function validateGeometry(obj, geoStruct)

    % geoStruct must be a struct
    assert(isstruct(geoStruct), 'Region of interest must be a struct.');
    
    % geoStruct must have a 'shape' field which is set to either 'rect' or 'ellipse'
    assert(isfield(geoStruct, 'shape'), 'geoStruct struct must have a ''shape'' field.');
    assert(ismember(geoStruct.shape, {'rect', 'ellipse'}),'geoStruct.shape can be either ''rect'' or ''ellipse''.');
    
    % geoStruct must have a 'units' field which is set to either 'degs' or 'microns'
    assert(isfield(geoStruct, 'units'), 'geoStruct struct must have a ''units'' field.');
    assert(ismember(geoStruct.units, {'degs', 'microns'}),'geoStruct.units can be either ''degs'' or ''microns''.');
    
    % geoStruct must have a 'center' field
    assert(isfield(geoStruct, 'center'), 'geoStruct struct must have a ''center'' field.');
    
    % geoStruct.center must have 2 elements (x,y)
    assert(numel(geoStruct.center) == 2, 'geoStruct.center must be a 2-element vector');
    
    if (strcmp(geoStruct.shape, 'rect'))
        % rect geoStruct must have 'width' and 'height' fields
        assert(isfield(geoStruct, 'width'), 'geoStruct struct must have a ''width'' field.');
        assert(isfield(geoStruct, 'height'), 'geoStruct struct must have a ''height'' field.');
    else
        % ellipse geoStruct must have 'rotation', 'minorAxisDiameter' and 'majorAxisDiameter' fields
        assert(isfield(geoStruct, 'rotation'), 'geoStruct struct must have a ''rotation'' field.');
        assert(isfield(geoStruct, 'minorAxisDiameter'), 'geoStruct struct must have a ''minorAxisDiameter'' field.');
        assert(isfield(geoStruct, 'majorAxisDiameter'), 'geoStruct struct must have a ''majorAxisDiameter'' field.');
    end
    
    % All OK, save the geometryStruct
    obj.geometryStruct = geoStruct;
end

