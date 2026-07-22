% s_coneEccentricities
%
% Estimate the cone apertures as a function of eccentricity.
% Uses the cMosaic
%
% This script generates plots to visualize the relationship between
% cone spacing, cone aperture, and cone aperture area with respect to
% eccentricity in degrees. The data is derived from a cMosaic object
% that simulates the spatial arrangement of cones in the retina.
%
% The script includes an option for rod intrusion adjustment.
%
% The source lattice follows the Watson eccentricity dependence but has
% approximately 10% smaller local cone spacing than the current Watson
% prediction. Its lattice generator uses a 2.15-micron minimum spacing,
% whereas the current Watson model predicts approximately 2.422 microns
% at the fovea. Regenerating the lattice with lambdaMinMicrons = 2.422
% would align the foveal spacing scale.
%
% See also
%   cMosaic

%% Whether to enable rod intrusion or not
ieInit;

rodIntrusion = true;

% Generate a narrow strip across the horizontal meridian. The narrow height
% retains local two-dimensional cone neighbors while avoiding the redundant
% off-meridian samples of the former 80-by-1-degree mosaic.
mosaicWidthDegs = 60;
mosaicHeightDegs = 0.25;
c = cMosaic(...
  'sourceLatticeSizeDegs', 64, ...
  'positionDegs', [0 0], ...
  'sizeDegs', [mosaicWidthDegs mosaicHeightDegs], ...
  'rodIntrusionAdjustedConeAperture', rodIntrusion, ...
  'coneApertureModifiers', struct('smoothLocalVariations', false), ...
  'useParfor', false);

coneEccentricity = sqrt(sum(c.coneRFpositionsDegs.^2,2)) .* sign(c.coneRFpositionsDegs(:,1));

% Evaluate the Watson cone-spacing model along the horizontal meridian.
% cMosaic imports the stored lattice after reversing its coordinates. Mirror
% the Watson x-coordinate here so its temporal and nasal predictions align
% with the displayed cMosaic coordinates; plot against the unmirrored
% watsonEccentricityDegs below.
watsonEccentricityDegs = linspace(-mosaicWidthDegs/2, mosaicWidthDegs/2, 1000)';
watsonEccentricityMicrons = -sign(watsonEccentricityDegs) .* 1000 .* ...
    RGCmodels.Watson.convert.rhoDegsToMMs(abs(watsonEccentricityDegs));
[watsonSpacingDegs, ~] = RGCmodels.Watson.compute.rfSpacingAtRetinalPositions( ...
    'right eye', [watsonEccentricityMicrons zeros(size(watsonEccentricityMicrons))], ...
    'cones', false);

%% Setup figure
ieFigure([], 'tall');
tiledlayout(3,1);

% Plot the variation in cone spacing with eccentricity
nexttile;
coneSpacingArcMin = c.coneRFspacingsDegs * 60; % Convert spacing to arc minutes
plot(coneEccentricity, coneSpacingArcMin, 'r.');
hold on;
plot(watsonEccentricityDegs, watsonSpacingDegs*60, 'k-', 'LineWidth', 1.5);
hold off;
set(gca, 'XTick', -30:5:30, 'XLim', mosaicWidthDegs/2 * [-1 1], 'YLim', [0 4], 'YTick', 0:0.5:10);
xtickangle(0);
grid on;
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
ylabel('cone spacing (arc min)');
legend({'cMosaic samples', 'Watson model'}, 'Location', 'northwest');
set(gca, 'FontSize', 20);

%% Plot the variation in cone aperture (microns) with eccentricity

nexttile
apertureMicrons = c.coneApertureDiametersDegs*c.micronsPerDegree;
plot(coneEccentricity, apertureMicrons, 'r.');
set(gca, 'XTick', -30:5:30, 'XLim', mosaicWidthDegs/2*[-1 1], 'YLim', [0 10], 'YTick', 0:2:10);
xtickangle(0)
grid on
ylabel('cone aperture (microns)');
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
set(gca, 'FontSize', 20);

%% Cone aperture area (um^2) assuming circular aperture (pi r2)

nexttile;
coneApertureArea = pi*(apertureMicrons/2).^2;
plot(coneEccentricity, coneApertureArea, 'r.');
set(gca, 'XTick', -30:10:30, 'XLim', mosaicWidthDegs/2*[-1 1], 'YLim', [0 50], 'YTick', 0:10:50);
xtickangle(0)
grid on
ylabel('cone aperture area (microns^2)');
xlabel('(  temporal  ) eccentricity (degs) (  nasal  )');
set(gca, 'FontSize', 24);

%%
