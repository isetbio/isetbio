% Create a 3D surface mesh with an image painted on it
%
% Use this to illustrate the type of output we might get from PBRT
% with the mesh implemented onto it.
%
% Demo for Joyce Liao about implications of curving the photoreceptor
% layer of the retina.  This might be from diseases in the retinal
% pigment epithelium that displace the photoreceptor layer.
%
% Amsler grid -  Used in diagnosis
% https://www.google.com/search?q=ansler+grid&rlz=1C5GCEM_enUS991US991&oq=ansler+grid&aqs=chrome..69i57j0i10i131i433i512l2j0i10i512l7.1889j0j4&sourceid=chrome&ie=UTF-8 
% 
% Delaunay triangulation
%
% TODO:  https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

%%  Build the 3D mesh
row = -50:50; col = -50:50;
[X, Y] = meshgrid(row,col);

% A function from X,Y to Z.  A circle for the moment.
%   X2 + Y2 + Z2 = radius2
radius = 500;
Z = -1*sqrt(radius^2 - X.^2 - Y.^2);
Z = Z - max(Z(:));

% Add a bump at some translated position T with some size sigma
sigma = 10;
T = [20 20];
gaussianBump = exp(-((X-T(1)).^2 + (Y-T(2)).^2)/(2*sigma^2));

sFactor = max(abs(Z(:)))*0.4;
Z = Z + sFactor*gaussianBump;
% visualize it
ieNewGraphWin;
mesh(X,Y,Z);

%%
img = imread('hatsC.jpg');
ieNewGraphWin; imagesc(img);

img2 = imresize(img,size(X));
ieNewGraphWin; imagesc(img2);

ieNewGraphWin;
s = mesh(X,Y,Z,img2);
s.FaceColor = 'flat';

%% Control the mesh in some interesting ways
sigma = 5;
T = [20 20];
gaussianBump = exp(-[(X-T(1)).^2 + (Y-T(2)).^2]/(2*sigma^2));
sFactor = 1e3;
s = mesh(X,Y,Z-sFactor*gaussianBump,img2);
s.FaceColor = 'flat';

%% Make some more test images, like the Amsler grid
hParams = harmonicP;
harmonic = imageHarmonic(hParams);
ieNewGraphWin; imagesc(harmonic)
harmonic = imresize(harmonic,size(X));
ieNewGraphWin;
s = mesh(X,Y,Z-sFactor*gaussianBump,harmonic);
s.FaceColor = 'flat';

%%
img = checkerboard(20,6,6);
img = imresize(img,size(X));
ieNewGraphWin; imagesc(img)
colormap(gray);

ieNewGraphWin;
s = mesh(X,Y,Z,img);
s.FaceColor = 'flat';

%%
scene = sceneCreate('letter');
letter = sceneGet(scene,'rgb');
ieNewGraphWin; imagesc(letter)
letter = imresize(letter,size(X));

ieNewGraphWin;
s = mesh(X,Y,Z-sFactor*gaussianBump,letter);
s.FaceColor = 'flat';

%%  This is not a good approach
% But we need to measure the distance along the surface of a mesh in
% Matlab.

% Here is a simple approach
%
% https://www.mathworks.com/support/search.html/answers/843605-geodesic-distances-on-a-curved-surface.html?fq%5B%5D=asset_type_name:answer&fq%5B%5D=category:matlab/computational-geometry&page=1 


% Two points in x,y.  
pt1 = [-3,2];
pt2 = [-1,1];

% Here is the surface
[x,y,z] = peaks(20);            % surface

% curve coordinates
x1 = linspace(pt1(1),pt1(2),20);
y1 = linspace(pt2(1),pt2(2),20);
z1 = interp2(x,y,z,x1,y1);

% Show the path.  It is not necessarily a geodesic, however.
surf(x,y,z,'edgecolor',[1 1 1]*0.8)
line(x1,y1,z1,'color','red','linewidth',2)

% The distance along the surface
L = sum(sqrt(diff(x1).^2+diff(y1).^2+diff(z1).^2))

% The distance along the plane underneath
D = norm(pt2 - pt1)
