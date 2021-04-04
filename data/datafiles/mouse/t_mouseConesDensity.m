%% t_mouseConesDensity
%
% Create some graphs comparing the density of the mouse cone photoreceptors
% with human.
%
% The mouse peak is about equal to the density at +/- 5 deg in human,
% around 16K or 18K cones/mm2. (% Number and Distribution of Mouse Retinal
% Cone Photoreceptors: Differences between an Albino (Swiss) and a
% Pigmented (C57/BL6) Strain)
%
% From the paper Geng et al. "Adaptive optics retinal imaging in the living
% mouse eye", every degree of visual angle corresponds to approximately 34
% µm for this model. In the human this is about 285 um.
%
% So for the equivalent stimulus (input referred) cone density of the mouse
% would need to have nearly 10x as high a cone density (linear).
%
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

%% Counting receptors
ieNewGraphWin;

mouse = 16*1e3*ones(size(coneDensity));

mid = ceil(size(X,2)/2);

plot(pos*1e3, coneDensity(:,mid),'b-','DisplayName','human','LineWidth',2);
hold on;
plot(pos*1e3, mouse(:,mid), 'k--','DisplayName','mouse','LineWidth',2);
plot(pos*1e3, mouse(:,mid)/10, 'r--', 'DisplayName','mouse corrected', 'LineWidth',2);
xlabel('Position (mm)')
ylabel('Density');
set(gca,'xscale','log');
grid on;
legend;

%%  Add a surface where the mouse density is

ieNewGraphWin;
[A, T] = meshgrid(ang, pos);
[X, Y] = pol2cart(A, T);

surf(X * 1e3, Y * 1e3, coneDensity); colormap(hsv)
xlabel('Position (mm)'); ylabel('Position (mm)'); 
zlabel('Cones / mm^2 ')

hold on;
mouseHigh = 16*1e3*ones(size(coneDensity))/10;   % Near the peak
mSurf = surf(X * 1e3, Y * 1e3, mouseHigh); 
mSurf.FaceAlpha = 0.3;
mSurf.FaceColor = [0.2 .2 .4];

%{
mouseLow = 6*1e3*ones(size(coneDensity))/10;     % Near the lowest
mSurf = surf(X * 1e3, Y * 1e3, mouseLow); 
mSurf.FaceAlpha = 0.3;
mSurf.FaceColor = [0.2 .4 .4];
%}

set(gca,'zscale','log');


%