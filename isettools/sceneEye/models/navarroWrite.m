function filename = navarroWrite(thisR)
% Write the navarro lens and support IoR files for a given accommodation
%
% Synopsis
%     filename = navarroWrite(thisR);
%
% Description
%   
%   The navarro_X_YY.dat file and associated index of refraction files
%   (iorX.spd) are written into the lens rendering directory. The navarro
%   model accounts for accommodation in the ior and lens files.  
% 
%   This is unclear to me.  What is the relationship between accommodation
%   and object distance?  I need to track this through the code.  But the
%   original text reads this way:
%
%   We use the 'object distance' slot to define the accommodation.  Until
%   the accommodation is less than 0.5 m, the impact of that factor is very
%   small on the IOR.
%
% Input
%  thisR:  The rendering recipe.  The accommodation (1 / focus distance)
%          will be specified in the lens file name.  The default was
%          navarro.dat. But maybe it should be navarro_X_YY_.dat
%
% Optional key/val pairs
%   N/A
%
% Outputs
%   filename:  Lens file
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
accom = (thisR.get('accommodation'));
na    = navarroLensCreate(accom);  % Diopters

% Build matrix and set focal Length
lensMatrix = [na.corneaA; na.corneaP; na.pupil; na.lensA; na.lensP];

% accom is in diopters (1/meters).  The base eye power is 60 diopters.  We
% add the accommodation to the base.
focalLength = 1 / (60.6061 + accom) * 10 ^ 3; % mm

%% Set up the lens sub-directory

lensDir = thisR.get('lens dir output');
if ~exist(lensDir,'dir'), mkdir(lensDir); end
lensFile = thisR.get('lens file');
% lensFile = fullfile(lensDir,'navarro.dat');

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

str = sprintf('\n# Accommodation (Diopters) %f \n',accom);
fprintf(fid, '%s', str);

str = '# See navarroLensCreate.m for adjusting accommodation';
fprintf(fid,'%s\n',str);
fclose(fid);

% The file should be there, so no warning should come from this.
thisR.set('lens file',lensFile);

%% Now write out the IoR files

% Our convention (which was hard coded in writeNavarroLensFile) is
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous
iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

% We assume the eye is accommodated to the object distance.  There is only
% a very very small impact of accommodation until the object is very close
% (less than 0.5 m).
[ior, wave]= navarroRefractiveIndices(accom);

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

fprintf('Wrote lens file to %s (accomm: %.2f D)\n',thisR.get('lensfile'),accom);

end
