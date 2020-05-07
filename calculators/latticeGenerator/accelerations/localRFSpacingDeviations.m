function spacingDeviations = localRFSpacingDeviations(rfPositions, tabulatedSpacing, tabulatedEcc)

    % Find distances to neighors
    neighborsNum = 1;
    spacings = localRFSpacings(rfPositions, neighborsNum);
    desiredSpacings = (lookUpValues(rfPositions, tabulatedEcc, tabulatedSpacing, neighborsNum))';
    spacingDeviations = (abs(desiredSpacings-spacings))./desiredSpacings;
end
