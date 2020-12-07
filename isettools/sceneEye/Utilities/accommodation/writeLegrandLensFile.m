function renderRecipe = writeLegrandLensFile(renderRecipe, workingFolder)
% DEPRECATED:  Write out the LeGrand lens file.
%
%    DEPRECATED:  Use the legrand* methods

% Syntax:
%   writeLegrandLensFile(renderRecipe, workingFolder)
%
% Description:
%    Write out the Le Grand theoreticl eye lens file. These values are
%    currently all fixed, but may be changeable in the future (I believe
%    the model has some simple accommodation models which we might
%    implement here.)
%
% Inputs:
%    renderRecipe  - ISET3d render @recipe
%    workingFolder - String. The working directory file path.
%
% Optional key/value pairs:
%    None.
%
% Outputs:
%    renderRecipe  - The modified Recipe
%
% See also
%  (BW fixed a lot of comments.  Started to add an example)

% History:
%   10/19/20  dhb Example checks for required iset3d function

% Examples:
%{
if (exist('piRecipeDefault','file'))
    thisR = piRecipeDefault;
    writeLegrandLensFile(thisR,thisR.get('output dir'));
else
    fprintf('This example requires iset3d on your path as well as isetbio\n');
end
%}

% Columns are: radiusX, radiusY, thickness, materialIndex, semiDiameter,
% conicConstantX, and conicConstantY
corneaA = [-7.8, -7.8, 0.55, 1, 4.820, 0, 0];
corneaP = [-6.5, -6.5, 3.05, 2, 4.341, 0, 0];
pupil = [0, 0, 0, 2, 2, 0, 0];
lensA = [-10.2, -10.2, 4, 3, 3.750, 0, 0];
lensP = [6, 6, 0, 4, 3.750, 0, 0];

%% Write out dispersion curves
% From Atchison, David A., and George Smith. "Chromatic dispersions of the
% ocular media of human eyes." JOSA A 22.1 (2005): 29-37.
wave = (400:10:800); % nm

% [Cornea; Aqueous; Lens; Vitreous]
n_inf = [1.3610 1.3221 1.3999 1.3208]';
K = [7.4147 7.0096 9.2492 6.9806]';
lambda_o = [130.0 130.0 130.0 130.0]';
V = [56 53 50 53]';

% Cornu dispersion equation
% Rows will be each ocular media
ior = n_inf + K ./ (wave - lambda_o);

iorNames = {sprintf('ior1_%0.2fdp_legrand.spd', 0), ...
    sprintf('ior2_%0.2fdp_legrand.spd', 0), ...
    sprintf('ior3_%0.2fdp_legrand.spd', 0), ...
    sprintf('ior4_%0.2fdp_legrand.spd', 0)};

% Write dispersion curves out to working folder
rtbWriteSpectrumFile(wave, ior(1,:), fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, ior(2,:), fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, ior(3,:), fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, ior(4,:), fullfile(workingFolder, iorNames{4}));

renderRecipe.camera.ior1.value = fullfile(workingFolder, iorNames{1});
renderRecipe.camera.ior2.value = fullfile(workingFolder, iorNames{2});
renderRecipe.camera.ior3.value = fullfile(workingFolder, iorNames{3});
renderRecipe.camera.ior4.value = fullfile(workingFolder, iorNames{4});

%% Write out lens file
lensFile = 'legrand.dat';
filename = fullfile(workingFolder, lensFile);
lensMatrix = [corneaA; corneaP; pupil; lensA; lensP];
focalLength = 1 / (59.94) * 10 ^ 3; % mm
fid = fopen(filename,'w');

str = sprintf('# Focal length (mm) \n');
fprintf(fid,'%s',str);
str = sprintf('%.3f\n',focalLength);
fprintf(fid,'%s',str);
str = sprintf(['# radiusX radiusY thickness materialIndex semiDiameter' ...
    ' conicConstantX conicConstantY\n']);
fprintf(fid,'%s',str);
for ii = 1:size(lensMatrix,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n', lensMatrix(ii,1), ...
        lensMatrix(ii, 2), lensMatrix(ii,3), lensMatrix(ii,4), ...
        lensMatrix(ii,5), lensMatrix(ii,6), lensMatrix(ii,7));
end
fclose(fid);

fprintf('Wrote out a new lens file: \n')
fprintf('%s \n \n', filename);

renderRecipe.camera.lensfile.value = filename;
renderRecipe.camera.lensfile.type = 'string';
        
end
