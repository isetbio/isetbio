function renderRecipe = writeLegrandLensFile(renderRecipe,workingFolder)
% Write out the Le Grand theoreticl eye lens file. These values are
% currently all fixed, but may be changeable in the future (I believe the
% model has some simple accommodation models which we might implement
% here.)
%
% Syntax:
%   writeLegrandLensFile(filename)
%
% Description:
%    Write out the Le Grand lens file
%
% Inputs:
%    filename - The filename to write to `
%
% Outputs:
%    None. 
%

% Columns are: radiusX, radiusY, thickness, materialIndex, semiDiameter,
% conicConstantX, and conicConstantY
corneaA = [-7.8  -7.8    0.55  1  4.820  0  0];
corneaP = [-6.5  -6.5    3.05  2  4.341  0  0];
pupil =   [ 0     0      0     2  2      0  0];
lensA =   [-10.2 -10.2   4     3  3.750  0  0];
lensP =   [ 6     6      0     4  3.750  0  0];

%% Write out dispersion curves

% The Le Grand eye doesn't seem to specify wavelength, so we'll use the
% same as the Arizona eye

wave = (400:10:800); % nm

% This is equivalent to IOR at 589.3 nm
% [cornea aqueous lens vitreous]
n_d = zeros(1,4);
n_d(1) = 1.377;
n_d(2) = 1.337;
n_d(3) = 1.42;
n_d(4) = 1.336;

% Abbe number
% [cornea aqueous lens vitreous]
V_d = zeros(1,4);
V_d(1) = 57.1;
V_d(2) = 61.3;
V_d(3) = 51.9;
V_d(4) = 61.1;

% Calculate dispersion curves
% n_f is IOR at 486.1 nm
% n_c is IOR at 656.3 nm
% V_d  = (n_d - 1)/(n_f - n_c)
% V_d = ( f(589.3) - 1 )/( f(486.1) - f(656.3) )
ior = cell(1,4); 

for ii = 1:length(n_d)
    
    m = (n_d(ii)-1)/(V_d(ii)*(486.1-656.3));
    b = n_d(ii)-m*589.3;
    ior{ii} = m*wave+b;        
end

iorNames = {sprintf('ior1_%0.2fdp_arizona.spd', 0), ...
    sprintf('ior2_%0.2fdp_arizona.spd', 0), ...
    sprintf('ior3_%0.2fdp_arizona.spd', 0), ...
    sprintf('ior4_%0.2fdp_arizona.spd', 0)};

% Write dispersion curves out to working folder
rtbWriteSpectrumFile(wave, ior{1}, fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, ior{2}, fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, ior{3}, fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, ior{4}, fullfile(workingFolder, iorNames{4}));

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
for ii=1:size(lensMatrix,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        lensMatrix(ii,1), lensMatrix(ii, 2), lensMatrix(ii,3), ...
        lensMatrix(ii,4), lensMatrix(ii,5), lensMatrix(ii,6), ...
        lensMatrix(ii,7));
end
fclose(fid);

fprintf('Wrote out a new lens file: \n')
fprintf('%s \n \n', filename);

renderRecipe.camera.lensfile.value = filename;
renderRecipe.camera.lensfile.type = 'string';
        
end
