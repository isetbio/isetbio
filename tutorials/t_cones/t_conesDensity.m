% Show how human cone density varies across the retina
%
% Description:
%    Plot the log10 cone density (cones / mm2) across the retinal surface.
%
%    The plot shows different angles and up to an eccentricity of 10 mm,
%    which is about +/-30 deg.
%

% History:
%    xx/xx/17  BW   ISETBIO Team, 2017
%    09/05/17  dhb  Clean up.
%    08/23/18  jnm  Formatting

%% Initialize
% ieInit;

%% Set range

pos = logspace(-5,-2,40);     % Spatial positions from fovea out (m)
ang = (0:.05:1.02) * 2 * pi;  % Angle around all 360 deg (plus a little)

% posDeg = logspace(-2, 1.5, 40);

%% Fill up a matrix with cone density

% This is cones per square millimeter
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

%% You could plot the individual lines
%{
	ieNewGraphWin;
	semilogx(pos * 1e3 * (1 / 3), coneDensity)
	xlabel('Position (deg)')
	ylabel('Density');
	grid on;
%}

%% But the surface plot looks snazzy
ieNewGraphWin;
[A, T] = meshgrid(ang, pos);
[X, Y] = pol2cart(A, T);

surf(X * 1e3, Y * 1e3, log10(coneDensity)); colormap(hsv)
xlabel('Position (mm)'); ylabel('Position (mm)'); 
zlabel('log_{10} Cones / mm^2 ')

%% Now show it with linear units

ieNewGraphWin;

surf(X * 1e3, Y * 1e3, coneDensity); colormap(hsv)
xlabel('Position (mm)'); ylabel('Position (mm)'); 
zlabel('Cones / mm^2 ')


%% END



