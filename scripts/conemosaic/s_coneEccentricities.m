% s_coneEccentricities
% 
% Obtain the cone eccentricities (signed) from a cmosaic
%

%% Whether to enable rod intrusion or not
rodIntrusionAdjustedConeAperture = true;

% Generate a 80x1 deg @cMosaic
c = cMosaic(...
  'sourceLatticeSizeDegs', 64, ...
  'positionDegs', [0 0], ...
  'sizeDegs', [80 1.0], ...
  'rodIntrusionAdjustedConeAperture', rodIntrusionAdjustedConeAperture);

theConeSignedEccentricityDegs = sqrt(sum(c.coneRFpositionsDegs.^2,2)) .* sign(c.coneRFpositionsDegs(:,1));

%% Setup figure
ieNewGraphWin([],'tall');
tiledlayout(3,1);

% Plot the variation in cone spacing with eccentricity
nexttile
plot(theConeSignedEccentricityDegs, c.coneRFspacingsDegs*60, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 4], 'YTick', 0:0.5:10);
xtickangle(0)
grid on
xlabel('(temporal retina) eccentricity (degs) (nasal retina)');
ylabel('cone spacings (arc min)');
set(gca, 'FontSize', 20);

%% Plot the variation in cone aperture with eccentricity
nexttile
plot(theConeSignedEccentricityDegs, c.coneApertureDiametersDegs*c.micronsPerDegree, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 10], 'YTick', 0:2:10);
xtickangle(0)
grid on
ylabel('cone apertures (arc min)');
xlabel('(temporal retina) eccentricity (degs) (nasal retina)');
set(gca, 'FontSize', 20);

%% Cone aperture area (um^2)

nexttile;
coneApertureArea = (c.coneApertureDiametersDegs*c.micronsPerDegree).^2;
plot(theConeSignedEccentricityDegs, coneApertureArea, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 50], 'YTick', 0:10:50);
xtickangle(0)
grid on 
ylabel('cone apertures (um^2)');
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
set(gca, 'FontSize', 20);


%% Export PDFs
%{
if (~rodIntrusionAdjustedConeAperture)
   NicePlot.exportFigToPDF('withoutRodIntrusion.pdf', hFig, 300);
else
   NicePlot.exportFigToPDF('withRodIntrusion.pdf', hFig, 300);
end
%}

%%