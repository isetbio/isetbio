function filename = arizonaWrite(thisR)
% Write out the Arizona lens file for a given accomodation
%
% Syntax:
%   arizonaWrite(thisR)
%
% Description:
%    Write out the Arizona lens file for a given accomodation
%
% Inputs:
%    thisR:   - rendering recipe
%
% Optional key/values:
%    None.
%
% Outputs:
%    filename:  lens file name
%
% See also
%   arizonaLensCreate; navarroWrite, legrandWrite
%
% Check with TL about the IOR files for the Arizona Eye Model
%

% History:
%   10/19/20  dhb Example checks for required iset3d function

% Examples:
%{
if (exist('piRecipeDefault','file'))
    myScene = sceneEye('slantedBar','human eye','arizona');
    thisR = myScene.get('recipe');
    arizonaWrite(thisR);
else
    fprintf('This example requires iset3d on your path as well as isetbio\n');
end
%}

%% Get the parameters given the accommodation

accommodation = thisR.get('accommodation');
az = arizonaLensCreate(accommodation);

%% Build matrix
lensMatrix = [az.corneaA; az.corneaP; az.pupil; az.lensA; az.lensP];

focalLength = 1 / (60.0 + accommodation) * 10 ^ 3; % mm

%% Set up the filename

lensDir = fullfile(thisR.get('output dir'),'lens');
if ~exist(lensDir,'dir'), mkdir(lensDir); end
filename = fullfile(lensDir,'arizona.dat');

%% Write the data into the file

fid = fopen(filename, 'w');

str = sprintf('# Focal length (mm) \n');
fprintf(fid, '%s', str);

str = sprintf('%.3f\n', focalLength);
fprintf(fid, '%s', str);

str = sprintf(['# radiusX radiusY thickness materialIndex semiDiameter' ...
    ' conicConstantX conicConstantY\n']);
fprintf(fid, '%s', str);

for ii = 1:size(lensMatrix, 1)
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        lensMatrix(ii, 1), lensMatrix(ii, 2), lensMatrix(ii, 3), ...
        lensMatrix(ii, 4), lensMatrix(ii, 5), lensMatrix(ii, 6), ...
        lensMatrix(ii, 7));
end

str = sprintf('\n# Accommodation (Diopters) %f \n',accommodation);
fprintf(fid, '%s', str);

fclose(fid);

%% Make sure the lens file is set properly

thisR.set('lens file',filename);


%% Now write out the IoR files for Arizona
%
% Waiting for TL comment

% Our convention (which was hard coded in writeNavarroLensFile) is
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous
iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

% We assume the eye is accommodated to the object distance.  There is only
% a very very small impact of accommodation until the object is very close
% (less than 0.5 m).
[ior,wave] = arizonaRefractiveIndices(accommodation);

% We will put these files next to the lens file (navarro.dat).
nSamples = numel(wave);
for ii=1:4
    filename = fullfile(thisR.get('lens dir output'),iorNames{ii});
    fid = fopen(filename, 'w');
    for jj = 1:nSamples
        fprintf(fid, '%d %f\n', wave(jj), ior(jj,ii));
    end
    fclose(fid);
    
    % Update the recipe with the ior files
    [~,str,~] = fileparts(filename);
    thisR.set(str,filename);
end

fprintf('Wrote lens file to %s (accomm: %.2f D)\n',thisR.get('lensfile'),accommodation);

end
