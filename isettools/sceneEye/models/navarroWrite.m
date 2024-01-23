function filename = navarroWrite(thisR,accommodation)
% Write the navarro lens and support IOR files for a given accommodation
%
% Synopsis
%     filename = navarroWrite(thisR,[accommodation]);
%
% Description
%   Write the navarro lens file and associated index of refraction (IOR)
%   files into the lens rendering directory. The navarro model accounts for
%   accommodation by changing the ior and lens files.
% 
% Input
%  thisR:  The rendering recipe.  The accommodation (1 / focus distance)
%          is specified in the lens file comment, and as of May 2023 it is
%          also stored in the focaldistance parameter. It can be returned
%          using
%
%               thisR.get('accommodation')
%
% Optional
%    accommodation - If not passed, we use thisR.get('accommodation')
%
% Outputs
%   filename:  lens file name.  
%       The file names include three digits that specify the accommodation
%       (diopters). So navarro-ABC.dat, and similarly for ior{1-4}-ABC.dat
%       and the accommodation is AB.C.
%
% See also
%   navarroLensCreate, navarroRefractiveIndices, arizonaWrite, legrandWrite

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

% TL believed that the usual accomodation value and the equivalent
% value for the Navarro model differ.  She inserted this function so
% that the Navarro model accommodation matched her Zemax calculation.
% See the comments in the function.
%
% The script t_eyeAccommodation images an edge, sweeping through
% accommodation from close, in focus, too far. The color fringe shifts
% as focus shift from too close to too far.
accommodation = convertToNavarroAccomm(accommodation);

% Round to 1 decimal place
accommodation = round(accommodation*10)/10;

na    = navarroLensCreate(accommodation);  % Diopters

% Build matrix and set focal Length
lensMatrix = [na.corneaA; na.corneaP; na.pupil; na.lensA; na.lensP];

% accommodation is in diopters (1/meters).  The base eye power is 60
% diopters.  We add the accommodation to the base.
focalLength = 1 / (60.6061 + accommodation) * 10 ^ 3; % mm

%% Set up the lens sub-directory

lensDir = thisR.get('lens dir output');
if ~exist(lensDir,'dir'), mkdir(lensDir); end
baseName = sprintf('navarro-%03.0f.dat',10*accommodation);

% Make the lens file for this accommodation and put it in the lens
% directory.
lensFile = fullfile(lensDir,baseName);
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
thisR.set('lens file',fullfile('lens',baseName));

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
iorNames = {'ior1','ior2','ior3','ior4'};

% We assume the eye is accommodated to the object distance.  There is only
% a very very small impact of accommodation until the object is very close
% (less than 0.5 m).
[ior, wave]= navarroRefractiveIndices(accommodation);

% We will put these files next to the lens file (navarro.dat).
nSamples = numel(wave);
for ii=1:4
    baseName = sprintf('%s-%03.0f.spd',iorNames{ii},accommodation);
    filename = fullfile(thisR.get('lens dir output'),baseName);
    fid = fopen(filename, 'w');
    for jj = 1:nSamples
        fprintf(fid, '%d %f\n', wave(jj), ior(jj,ii));
    end
    fclose(fid);
    
    % Update the recipe with the ior files
    filename = fullfile('lens',baseName);
    thisR.set(iorNames{ii},filename);
end

fprintf('Wrote lens file to %s (navarro adjusted accomm: %.2f D)\n',thisR.get('lensfile'),accommodation);

end
