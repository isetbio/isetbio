function [] = writeLegrandLensFile(filename)
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
corneaA = [-7.8  -7.8    0.55  1  4.820  -0.26  -0.26];
corneaP = [-6.5  -6.5    3.05  2  4.341  -0.26  -0.26];
pupil =   [ 0     0      0     2  2       0      0   ];
lensA =   [-10.2 -10.2   4     3  3.750  -3.136 -3.136];
lensP =   [ 6     6      0     4  3.750  -1     -1];

%% Write out dispersion curves

% TODO: Where did I get these from? I think it might have been an Atchinson
% paper somewhere...
% Columns: [cornea aqueous lens vitreous]
wave = 400:10:800;
ior = ieReadSpectra(fullfile(piRootPath,'data','lens','IORofEye.mat'),wave);

iorNames = {sprintf('ior1_%0.2fdp.spd', 0), ...
    sprintf('ior2_%0.2fdp.spd', 0), ...
    sprintf('ior3_%0.2fdp.spd', 0), ...
    sprintf('ior4_%0.2fdp.spd', 0)};

% Get working directory from filename
[workingFolder,~,~] = fileparts(filename);

% Write dispersion curves out to working folder
rtbWriteSpectrumFile(wave, ior(:,1), fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, ior(:,2), fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, ior(:,3), fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, ior(:,4), fullfile(workingFolder, iorNames{4}));

%% Build matrix
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

end
