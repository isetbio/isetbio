function regenerateConePositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile, varargin)

    p = inputParser;
    p.addParameter('customDegsToMMsConversionFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customMMsToDegsConversionFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customRFspacingFunction', [], @(x) (isempty(x) || isa(x,'function_handle')));
    p.addParameter('customMinRFspacing', [], @(x) (isempty(x) || isscalar(x)));
    p.parse(varargin{:});
    
    % Regenerate lattice whose FOV is large enough to encopass the desired size at the desired eccentricity
    fovDegs = (sqrt(sum(obj.eccentricityDegs.^2,2)) + max(obj.sizeDegs))*1.3;

    obj.coneRFpositionsMicrons = retinalattice.generatePatch(fovDegs, ...
        'cones', obj.whichEye, exportHistoryToFile, visualizeConvergence, obj.useParfor, maxIterations, ...
        'randomSeed', obj.randomSeed, ...
        'customDegsToMMsConversionFunction', p.Results.customDegsToMMsConversionFunction, ...
        'customRFspacingFunction', p.Results.customRFspacingFunction, ...
        'customMinRFspacing', p.Results.customMinRFspacing);
    
    % Convert to degs
    if (~isempty(p.Results.customMMsToDegsConversionFunction))
        obj.coneRFpositionsDegs = p.Results.customMMsToDegsConversionFunction(obj.coneRFpositionsMicrons*1e-3);
    else
        obj.coneRFpositionsDegs = RGCmodels.Watson.convert.rhoMMsToDegs(obj.coneRFpositionsMicrons*1e-3);
    end
    
    % Compute spacings (which determine apertures)
    obj.coneRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsDegs);
    obj.coneRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.coneRFpositionsMicrons);
    
    % Crop to desired ROI in degs
    diff = abs(bsxfun(@minus, obj.coneRFpositionsDegs, obj.eccentricityDegs));
    idx = find((diff(:,1) <= 0.5*obj.sizeDegs(1)) & (diff(:,2) <= 0.5*obj.sizeDegs(2)));
    obj.coneRFpositionsDegs = obj.coneRFpositionsDegs(idx,:);
    obj.coneRFpositionsMicrons = obj.coneRFpositionsMicrons(idx,:);
    obj.coneRFspacingsDegs = obj.coneRFspacingsDegs(idx);
    obj.coneRFspacingsMicrons = obj.coneRFspacingsMicrons(idx);
end
