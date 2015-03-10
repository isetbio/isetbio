% s_XYZilluminantTransforms
%
%  This script illustrates a method for calculating diagonal and 3x3
%  transforms that map the XYZ values for surfaces under on light into the
%  XYZ values under another light
% 
%  Create a scene that consists of patches illuminated by light A
%       these patches can be created by randomly sampling surfaces reflectances
%       the selection of surfaces will matter
%  Create another scene that consists of the same patches illuminated by light B
%  Calculate a 3x3 matrix that maps the XYZ values for surfaces under light A
%       into the XYZ values of surfaces under light B
%  Calculate a diagonal matrix that maps the XYZ values for surfaces under light A
%       into the XYZ values of surfaces under light B
%  Notice that the 3x3 does a better job than the diagonal (as expected)
% 
% Copyright ImagEval Consultants, LLC, 2012.

%% Randomly select reflectances
% The files containing the reflectances are in ISET format, readable by 
% s = ieReadSpectra(sFiles{1});
sFiles = cell(1,4);
sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','Food_Vhrel.mat');
sFiles{3} = fullfile(isetRootPath,'data','surfaces','reflectances','DupontPaintChip_Vhrel.mat');
sFiles{4} = fullfile(isetRootPath,'data','surfaces','reflectances','HyspexSkinReflectance.mat');

% The number of samples from each of the data sets, respectively
sSamples = [12,12,24,24];    % 

% How many row/col spatial samples in each patch (they are square)
pSize = 24;    % Patch size
wave =[];      % Whatever is in the file
grayFlag = 0;  % No gray strip
sampling = 'no replacement';

scene = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayFlag,sampling);

% Show it on the screen
vcAddAndSelectObject(scene); sceneWindow;

%% Change the illumination
% Change from the default illuminant (equal energy) to D65

wave = sceneGet(scene,'wave');  d65 = ieReadSpectra('D65',wave);
sceneD65 = sceneAdjustIlluminant(scene,d65);
sceneD65 = sceneSet(sceneD65,'name','Reflectance Chart D65');
vcAddAndSelectObject(sceneD65); sceneWindow;

%% Add a gray strip column
grayStrip = 1;
sceneGray = sceneReflectanceChart(sFiles,sSamples,pSize,wave,grayStrip);
sceneGray = sceneSet(sceneGray,'name','Reflectance Chart EE Gray Strip');
vcAddAndSelectObject(sceneGray); sceneWindow;

%%  Solve for matrix relating the chart under two different lights

% This are the surfaces under a D65 light
xyz1 = sceneGet(scene,'xyz');

% This is a nSample x 3 representation of the surfaces under D65
xyz1 = RGB2XWFormat(xyz1);

% This are the surfaces under a Tungsten light
scene2 = sceneAdjustIlluminant(scene,'Tungsten.mat');
xyz2 = sceneGet(scene2,'xyz');

% This is a nSample x 3 representation of the surfaces under Tungsten
xyz2 = RGB2XWFormat(xyz2);

% We are looking for a 3x3 matrix, L, that maps
%
%    xyz1 = xyz2 * L
%    L = inv(xyz2'*xyz2)*xyz2'*xyz1 = pinv(xyz2)*xyz1
%
% Or, we just use the \ operator from Matlab for which inv(A)*B is A\B
L = xyz2 \ xyz1;

% To solve with just a diagonal, do it one column at a time
D = zeros(3,3);
for ii=1:3
    D(ii,ii) = xyz2(:,ii) \ xyz1(:,ii);
end

%% Plot predicted versus actual
% vcNewGraphWin; pred2 = xyz2*L; plot(xyz1(:),pred2(:),'.')
% vcNewGraphWin; pred2 = xyz2*D; plot(xyz1(:),pred2(:),'.')


%% End

