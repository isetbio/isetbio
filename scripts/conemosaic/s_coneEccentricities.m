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
hFig = figure(2); clf;
set(hFig, 'Position', [10 10 950 1250], 'Color', [1 1 1]);

% Plot the variation in cone spacing with eccentricity
subplot(2,1,1)
plot(theConeSignedEccentricityDegs, c.coneRFspacingsDegs*60, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 4], 'YTick', 0:0.5:10);
xtickangle(0)
grid on
xlabel('(temporal retina) eccentricity (degs) (nasal retina)');
ylabel('cone spacings (arc min)');
set(gca, 'FontSize', 20);

%% Plot the variation in cone aperture with eccentricity
subplot(2,1,2)
plot(theConeSignedEccentricityDegs, c.coneApertureDiametersDegs*60, 'r.');
set(gca, 'XTick', -40:5:40, 'XLim', 40*[-1 1], 'YLim', [0 4], 'YTick', 0:0.5:10);
xtickangle(0)
grid on
ylabel('cone apertures (arc min)');
xlabel('(temporal retina) eccentricity (degs) (nasal retina)');
set(gca, 'FontSize', 20);

% Export PDFs
if (~rodIntrusionAdjustedConeAperture)
   NicePlot.exportFigToPDF('withoutRodIntrusion.pdf', hFig, 300);
else
   NicePlot.exportFigToPDF('withRodIntrusion.pdf', hFig, 300);
end