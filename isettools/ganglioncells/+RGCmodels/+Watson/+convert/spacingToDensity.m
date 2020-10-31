function density = spacingToDensity(spacing)
% Convert spacing to density in a perfect hex mosaic (Equation A4)
    density = 2.0./(sqrt(3.0)*spacing.^2);
end