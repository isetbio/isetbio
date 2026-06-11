function theSmoothedSpacings = smoothSpacings(rfSpacings, nearbyRFindices)

    % Smooth variations in RF spacings
    theSmoothedSpacings = 0*rfSpacings;
    
    parfor rfIndex = 1:numel(rfSpacings)
        theSmoothedSpacings(rfIndex) = median(rfSpacings(nearbyRFindices(rfIndex,:)));
    end
end
