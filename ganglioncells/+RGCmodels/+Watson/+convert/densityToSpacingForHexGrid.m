function spacing = densityToSpacingForHexGrid(density)
% Convert density to spacing in a perfect hex mosaic (Equation A4)
    spacing = sqrt(2.0./(sqrt(3.0)*density));
end

