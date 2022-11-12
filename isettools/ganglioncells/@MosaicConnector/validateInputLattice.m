function validateInputLattice(obj, theLattice, latticeName)

    % Must be a struct
    assert(isstruct(theLattice), 'The %s lattice must be a struct.', latticeName);
    
    % Must have field
    assert(isfield(theLattice, 'name')&&(ischar(theLattice.name)), 'The %s lattice must have a ''name'' field.', latticeName);
    assert(isfield(theLattice, 'DegsToMMsConversionFunction')&&(isa(theLattice.DegsToMMsConversionFunction,'function_handle')), 'The %s lattice must have a ''DegsToMMsConversionFunction'' function handle field.', latticeName);
    assert(isfield(theLattice, 'MMsToDegsConversionFunction')&&(isa(theLattice.MMsToDegsConversionFunction,'function_handle')), 'The %s lattice must have a ''MMsToDegsConversionFunction'' function handle field.', latticeName);
    assert(isfield(theLattice, 'RFpositionsMicrons'), 'The %s lattice must have an ''RFpositionsMicrons'' field.', latticeName);
    
    % Compute spacings
    [theLattice.RFspacingsMicrons, nearbyRFindices] = RGCmodels.Watson.convert.positionsToSpacings(theLattice.RFpositionsMicrons);
    theLattice.nearbyRFindices = nearbyRFindices';

    assert(isfield(theLattice, 'RFspacingsMicrons'), 'The %s lattice must have an ''RFspacingsMicrons'' field.', latticeName);
    % Matrix dimensions must be valid
    [nPos,dimensions] = size(theLattice.RFpositionsMicrons);
    
    assert(dimensions == 2, ...
        'The ''RFpositionsMicrons'' field of the %s lattice must be an N x 2 matrix. The passed data is %d x %d', latticeName, nPos, dimensions);
    [nPos2,dimensions2] = size(theLattice.RFspacingsMicrons);
    if (nPos2 == 1)
        theLattice.RFspacingsMicrons = theLattice.RFspacingsMicrons(:);
        [nPos2,dimensions2] = size(theLattice.RFspacingsMicrons);
    end
    
    assert(dimensions2 == 1, ...
        'The ''RFspacingsMicrons'' field of the %s lattice must be an N x 1 matrix. The passed data is %d x %d', latticeName, nPos2, dimensions2);
    assert(nPos == nPos2, 'The ''RFpositionsMicrons'' field of the %s lattice does not have the same rows (%d) as the ''RFspacingsMicrons'' field (%d)', latticeName,nPos, nPos2);
    
    switch (latticeName)
        case 'source'
            % Smooth source lattice spacings
            if (obj.smoothSourceLatticeSpacings)
                theLattice.RFspacingsMicrons = MosaicConnector.smoothSpacings(...
                    theLattice.RFspacingsMicrons, theLattice.nearbyRFindices);
            end
            obj.sourceLattice = theLattice;

        case 'destination'
            % Smooth destination lattice spacings
            if (obj.smoothDestinationLatticeSpacings)
                theLattice.RFspacingsMicrons = MosaicConnector.smoothSpacings(...
                    theLattice.RFspacingsMicrons, theLattice.nearbyRFindices);
            end
            obj.destinationLattice = theLattice;

        otherwise
            error('lattice name (''%s'') is invalid. It must be either ''source'' or ''destination''.', latticeName);
    end

    

end