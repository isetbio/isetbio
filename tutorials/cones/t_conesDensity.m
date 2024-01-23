% Human cone density varies across the retina
%
% Description:
%    Plot the cone density across the retinal surface.
%
%    The Curcio plot shows different angles and up to an eccentricity
%    of 10 mm, which is about +/-30 deg.  The Song plots are up to 2
%    mm.
%
%    Sources for cone density estimate
%        'Curcio1990'         From Figure 6 of Ref 1 below (default).
%        'Song2011Old'        From Table 1 of Ref 2 below, old subjects data.
%        'Song2011Young'      From Table 1 of Ref 2 below, young subjects data.
%
%   1) Curcio, C. A., Sloan, K. R., Kalina, R. E. and Hendrickson, A. E.
%      (1990), Human photoreceptor topography. J. Comp. Neurol., 292:
%      497?523. doi: 10.1002/cne.902920402
%   2) Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., & Burns, S. A.
%      (2011). Variation of Cone Photoreceptor Packing Density with Retinal
%      Eccentricity and Age. Investigative Ophthalmology & Visual Science,
%      52(10), 7376-7384. http://doi.org/10.1167/iovs.11-7199
%
% See also
%   t_conesPhotopigment, t_lensTransmittance
%

% History:
%    xx/xx/17  BW   ISETBIO Team, 2017
%    09/05/17  dhb  Clean up.
%    08/23/18  jnm  Formatting

%% Initialize
ieInit;

%% Set range

pos = logspace(-5,-2,40);     % Spatial positions from fovea out (mm)
ang = (0:.1:1.02) * 2 * pi;  % Angle around all 360 deg (plus a little)

%% Fill up a matrix with cone density

% Possible sources are:
%    'Curcio1990' 'Song2011Old' 'Song2011Young'

coneDensity = zeros(length(pos), length(ang));
for jj = 1:length(ang)
    for ii = 1:length(pos)
        coneDensity(ii, jj) = coneDensityReadData(...
            'eccentricity', pos(ii), 'eccentricityUnits', 'm',...
            'angle', ang(jj), 'angleUnits','deg', ...
            'coneDensitySource','Curcio1990',...
            'whichEye', 'left');
    end
end

%% You can plot an individual line of density vs. position

ieNewGraphWin;
semilogx(pos * 1e3 * (1 / 3), coneDensity)
xlabel('Position (deg)')
ylabel('Density');
grid on;

%% But the surface plot looks snazzy
ieNewGraphWin;
[A, T] = meshgrid(ang, pos);
[X, Y] = pol2cart(A, T);
surf(X * 1e3, Y * 1e3, log10(coneDensity)); colormap(hsv)
xlabel('Position (mm)'); ylabel('Position (mm)'); 
zlabel('log_{10} Cones / mm^2 ')

%% Plotted in linear units makes the significance of the fovea very clear

surf(X * 1e3, Y * 1e3, coneDensity); colormap(hsv)
xlabel('Position (mm)'); ylabel('Position (mm)'); 
zlabel('Cones / mm^2 (log)')
set(gca,'zscale','log');
%% END



