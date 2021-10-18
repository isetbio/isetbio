% Static method to validate an ROI struct
function validateROI(roi)
    % roi must be a struct
    assert(isstruct(roi), 'Region of interest must be a struct.');
    
    % roi must have a 'shape' field which is set to either 'rect' or 'ellipse'
    assert(isfield(roi, 'shape'), 'roi struct must have a ''shape'' field.');
    assert(ismember(roi.shape, {'rect', 'ellipse'}),'roi.shape can be either ''rect'' or ''ellipse''.');
    
    % roi must have a 'units' field which is set to either 'degs' or 'microns'
    assert(isfield(roi, 'units'), 'roi struct must have a ''units'' field.');
    assert(ismember(roi.units, {'degs', 'microns'}),'roi.units can be either ''degs'' or ''microns''.');
    
    % roi must have a 'center' field
    assert(isfield(roi, 'center'), 'roi struct must have a ''center'' field.');
    
    % roi.center must have 2 elements (x,y)
    assert(numel(roi.center) == 2, 'roi.center must be a 2-element vector');
    
    if (strcmp(roi.shape, 'rect'))
        % rect roi must have 'width' and 'height' fields
        assert(isfield(roi, 'width'), 'roi struct must have a ''width'' field.');
        assert(isfield(roi, 'height'), 'roi struct must have a ''height'' field.');
    else
        % ellipse roi must have 'rotation', 'minorAxisDiameter' and 'majorAxisDiameter' fields
        assert(isfield(roi, 'rotation'), 'roi struct must have a ''rotation'' field.');
        assert(isfield(roi, 'minorAxisDiameter'), 'roi struct must have a ''minorAxisDiameter'' field.');
        assert(isfield(roi, 'majorAxisDiameter'), 'roi struct must have a ''majorAxisDiameter'' field.');
    end
end