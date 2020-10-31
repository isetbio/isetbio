function thisR = legrandWrite(thisR)
% Write out the Le Grand lens file and support IoR files in the output directory
%
% Syntax:
%   legrandWrite(thisR)
%
% Description:
%    Write out the Le Grand theoretical eye lens file. These values are
%    currently all fixed, but may be changeable in the future (I believe
%    the model has some simple accommodation models which we might
%    implement here.)
%
% Inputs:
%    thisR   - Render recipe
%
% Optional key/value pairs:
%    None.
%
% Outputs:
%    thisR  - The updated Recipe
%
% See also
%   navarroWrite

% History:
%   10/19/20  dhb Example checks for required iset3d function

% Examples:
%{
if (exist('piRecipeDefault','file'))
    thisR = piRecipeDefault;
    thisR.camera = piCameraCreate('humaneye','lens file','legrand.dat');
    legrandWrite(thisR);
else
    fprintf('This example requires iset3d on your path as well as isetbio\n');
end
%}
%% Parameters

% This returns the model parameters for the geometric components of the
% LeGrand model eye.
lg = legrandLensCreate;

% Assemble them
lensMatrix = [lg.corneaA; lg.corneaP; lg.pupil; lg.lensA; lg.lensP];

% Express the focal length
focalLength = 1 / (59.94) * 10 ^ 3; % mm

%% Write out lens file into the lens subdirectory.

% Check that the directory is there
lensDir = fullfile(thisR.get('output dir'),'lens');
if ~exist(lensDir,'dir')
    fprintf('Making lens directory\n');
    mkdir(lensDir); 
end

% Set up the file name
filename = fullfile(lensDir,'legrand.dat');

% Open and write
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

%% The camera should be a realistic eye type

c = thisR.get('camera subtype');
if ~isequal(c, 'realisticEye')
    warning('LeGrand lens written, but camera subtype is %s\n',c);
end

% Store the name
thisR.set('lens file',filename);

%% Write out the IOR files in the same directory

[ior,wave] = legrandRefractiveIndices;
iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

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
        
fprintf('Wrote lens file to %s\n',thisR.get('lensfile'));

end
