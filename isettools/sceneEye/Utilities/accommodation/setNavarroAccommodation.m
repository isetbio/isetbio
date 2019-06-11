function renderRecipe = setNavarroAccommodation(...
    renderRecipe, accommodation, workingFolder)
% Change renderRecipe to match the accommodation
%
% Syntax:
%   renderRecipe = setNavarroAccommodation(..
%       renderRecipe, accommodation, workingFolder)
%
% Description:
%    We change the fields of the renderRecipe to match accommodation. As
%    accommodation changes, the lens file will change, as will the index of
%    refraction for the lens media. We write these new files out and
%    reference them in the structure.
%
% Inputs:
%    renderRecipe  - Object. Un-modified renderRecipe object.
%    accommodation - Numeric. The accommodation to shape the modified
%                    renderRecipe by.
%    workingFolder - String. The file location to write the new
%                    renderRecipe to.
%
% Outputs:
%    renderRecipe  - Object. The modified renderRecipe.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/17  TL   Created by Trisha Lian IESTBIO Team 2017
%    12/19/17  jnm  Formatting
%    01/25/19  JNM  Changed output of filenames to remove extra periods,
%                   which is breaking Windows executions.
%    05/29/19  JNM  Second documentation pass (minor tweaks)

%% Check and make sure this recipe includes a realisticEye
if(~strcmp(renderRecipe.camera.subtype, 'realisticEye'))
    warning('The camera type is not a realisticEye. Returning untouched.');
    return;
end

%% Check inputs
if(~(accommodation >= 0 && accommodation <= 10))
    error('Accommodation must be between 0 and 10 diopters.');
end

if ~exist(workingFolder, 'dir')
    error('Working folder does not exist.');
end

%% Convert accommodation
% See the function description for more information on why this is needed.
navarroAccom = convertToNavarroAccomm(accommodation);

%% Write out ocular media spectra files
% We calculate the dispersion curves of each ocular media. In the lens
% file, each surface material is linked to an "ior slot" (ior1, 
% ior2, etc.) When the ray is traveling through that material, it will
% follow the curve defined by the spectrum in the corresponding slot.

% Our convention (hard coded in writeNavarroLensFile) is always:
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous
wave = (400:10:800); % nm
[cor, aqu, len, vit] = getNavarroRefractiveIndices(wave, accommodation);

accStr = sprintf('%0.2f', accommodation);
accStr = strrep(accStr, '.', '_');
iorNames = {sprintf('ior1_%sdp.spd', accStr), ...
    sprintf('ior2_%sdp.spd', accStr), ...
    sprintf('ior3_%sdp.spd', accStr), ...
    sprintf('ior4_%sdp.spd', accStr)};

rtbWriteSpectrumFile(wave, cor, fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, aqu, fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, len, fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, vit, fullfile(workingFolder, iorNames{4}));

renderRecipe.camera.ior1.value = fullfile(workingFolder, iorNames{1});
renderRecipe.camera.ior2.value = fullfile(workingFolder, iorNames{2});
renderRecipe.camera.ior3.value = fullfile(workingFolder, iorNames{3});
renderRecipe.camera.ior4.value = fullfile(workingFolder, iorNames{4});

%% Attach lens file and set retina radius
% For navarro, the lens file will change depending on accomodation. Here we
% can write it out to a file to be read in later.
lensFile = sprintf('navarroAccomodated_%s.dat', accStr);
writeNavarroLensFile(navarroAccom, fullfile(workingFolder, lensFile));
fprintf('Wrote out a new lens file: \n')
fprintf('%s \n \n', fullfile(workingFolder, lensFile));

renderRecipe.camera.lensfile.value = fullfile(workingFolder, lensFile);
renderRecipe.camera.lensfile.type = 'string';

end