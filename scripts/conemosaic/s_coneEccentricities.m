% s_coneEccentricities
%
% Generate plots to visualize the relationship between cone spacing,
% cone aperture, and cone aperture area with respect to eccentricity
% in degrees. 
% 
% The estimates are derived from a cMosaic object that simulates
% the spatial arrangement of cones in the retina.  The aperture
% diameters are calculated from the cone spacing, multiplying by a
% reduction based on these factors and potentially a rod intrusion
% factor.
%
% eccVaryingApertures = obj.coneRFspacingsMicrons * ...
%                      obj.coneDiameterToSpacingRatio * ...
%                      obj.coneApertureToDiameterRatio;
%
%
% See also
%   cMosaic, NicePlot.exportFigToPDF

%% Whether to enable rod intrusion or not
rodIntrusion = true;

% Generate a 80x1 deg @cMosaic
c = cMosaic(...
  'sourceLatticeSizeDegs', 64, ...
  'positionDegs', [0 0], ...
  'sizeDegs', [80 1.0], ...
  'rodIntrusionAdjustedConeAperture', rodIntrusion);

coneEccentricity = sqrt(sum(c.coneRFpositionsDegs.^2,2)) .* sign(c.coneRFpositionsDegs(:,1));

%% Setup figure
ieFigure([],'tall');
tiledlayout(3,1);

% Plot the variation in cone spacing with eccentricity
nexttile;
coneSpacingArcMin = c.coneRFspacingsDegs * 60; % Convert spacing to arc minutes
plot(coneEccentricity, coneSpacingArcMin, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40 * [-1 1], 'YLim', [0 4], 'YTick', 0:0.5:10);
xtickangle(0);
grid on;
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
ylabel('cone spacing (arc min)');
set(gca, 'FontSize', 20);

%% Plot the variation in cone aperture (microns) with eccentricity

nexttile
apertureMicrons = c.coneApertureDiametersDegs*c.micronsPerDegree;
plot(coneEccentricity, apertureMicrons, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 10], 'YTick', 0:2:10);
xtickangle(0)
grid on
ylabel('cone aperture (arc min)');
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
set(gca, 'FontSize', 20);

%% Cone aperture area (um^2) assuming circular aperture (pi r2)

nexttile;
coneApertureArea = pi*(apertureMicrons/2).^2;
p = plot(coneEccentricity, coneApertureArea, 'r.');
set(gca, 'XTick', -40:10:40, 'XLim', 40*[-1 1], 'YLim', [0 50], 'YTick', 0:10:50);
xtickangle(0)
grid on
ylabel('cone aperture (um^2)');
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
set(gca, 'FontSize', 24);


%% Export PDFs
%{
if (~rodIntrusionAdjustedConeAperture)
   NicePlot.exportFigToPDF('withoutRodIntrusion.pdf', hFig, 300);
else
   NicePlot.exportFigToPDF('withRodIntrusion.pdf', hFig, 300);
end
%}

%%