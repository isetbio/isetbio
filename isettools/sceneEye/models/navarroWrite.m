function filename = navarroWrite(thisR,accommodation)
% Write the navarro lens and support IoR files for a given accommodation
%
% Synopsis
%     filename = navarroWrite(thisR,[accommodation]);
%
% Description
%   Write the navarro.dat file and associated index of refraction
%   files (iorX.spd) into the lens rendering directory. The navarro
%   model accounts for accommodation by changing the ior and lens
%   files.  The retinalDistance does not change.
% 
% Input
%  thisR:  The rendering recipe.  The accommodation (1 / focus distance)
%          is specified in the lens file comment.  The file can be
%          read and the value extracted using
%
%                thisR.get('accommodation')
%
%  accom:  Specify an accommodation value.  Default is 0 (Accommodated
%  to infinity).  
%
% Outputs
%   filename:  Lens file, full path to navarro.dat
%
% See also
%   navarroLensCreate, navarroRefractiveIndices

% History:
%   10/19/20  dhb Example checks for required iset3d function

% Examples:
%{
if (exist('piRecipeDefault','file'))
    myScene = sceneEye('slantedBar');
    thisR = myScene.recipe;
    navarroWrite(thisR);
else
    fprintf('This example requires iset3d on your path as well as isetbio\n');
end
%}


%% Writes out the navarro.dat file in the lens directory of the output
if notDefined('accommodation')
    accommodation = (thisR.get('accommodation'));
end
na    = navarroLensCreate(accommodation);  % Diopters

% Build matrix and set focal Length
lensMatrix = [na.corneaA; na.corneaP; na.pupil; na.lensA; na.lensP];

% accommodation is in diopters (1/meters).  The base eye power is 60
% diopters.  We add the accommodation to the base.
focalLength = 1 / (60.6061 + accommodation) * 10 ^ 3; % mm

%% Set up the lens sub-directory

lensDir = thisR.get('lens dir output');
if ~exist(lensDir,'dir'), mkdir(lensDir); end

% Make the lens file for this accommodation and put it in the lens
% directory.
lensFile = fullfile(lensDir,'navarro.dat');
thisR.set('lens file',lensFile);

%% Do the writing
fid = fopen(lensFile, 'w');

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

str = '# See navarroLensCreate.m for adjusting accommodation';
fprintf(fid,'%s\n',str);
fclose(fid);

% The file should be there, so no warning should come from this.
thisR.set('lens file',lensFile);

%% Now write out the IoR files

% We calculate and write out the index of refraction curves of each ocular
% media. Each surface boundary is linked to an "ior slot" (ior1, ior2,
% etc.) When the ray is traveling through that material, it will follow the
% curve defined by the spectrum in the corresponding interface.
%
% Our convention (hard coded in writeNavarroLensFile) is always these
% interfaces: 
%
%   ior1 --> air-cornea
%   ior2 --> cornea-aqueuous
%   ior3 --> aqueous-lens
%   ior4 --> lens-vitreous
%
iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

% We assume the eye is accommodated to the object distance.  There is only
% a very very small impact of accommodation until the object is very close
% (less than 0.5 m).
[ior, wave]= navarroRefractiveIndices(accommodation);

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
