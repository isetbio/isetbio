function renderRecipe = setAccommodation(renderRecipe,accommodation,workingFolder)
% We change the fields of the renderRecipe to match accommodation. As
% accommodation changes, the lens file will change, as will the index of
% refraction for the lens media. We write these new files out and reference
% them in the structure.

% Trisha Lian IESTBIO Team 2017

%% Check and make sure this recipe includes a realisticEye

if(~strcmp(renderRecipe.camera.subtype,'realisticEye'))
    warning('The camera type is not a realisticEye. Returning untouched.');
    return;
end

%% Check inputs
if(~(accommodation >= 0 && accommodation <= 10))
    error('Accommodation must be between 0 and 10 diopters.');
end

if(~exist(workingFolder,'dir'))
    error('Working folder does not exist.');
end

%% Convert accommodation
% See the description of this function for more information on why this is
% necessary. 
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

wave = (400:10:800); % um

[cor,aqu,len,vit] = getNavarroRefractiveIndices(wave,accommodation);

iorNames = {sprintf('ior1_%0.2fdp.spd',accommodation),...
    sprintf('ior2_%0.2fdp.spd',accommodation),...
    sprintf('ior3_%0.2fdp.spd',accommodation),...
    sprintf('ior4_%0.2fdp.spd',accommodation)};

rtbWriteSpectrumFile(wave,cor,fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave,aqu,fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave,len,fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave,vit,fullfile(workingFolder, iorNames{4}));

renderRecipe.camera.ior1.value = fullfile(workingFolder,iorNames{1});
renderRecipe.camera.ior2.value = fullfile(workingFolder,iorNames{2});
renderRecipe.camera.ior3.value = fullfile(workingFolder,iorNames{3});
renderRecipe.camera.ior4.value = fullfile(workingFolder,iorNames{4});

%% Attach lens file and set retina radius

% For navarro, the lens file will change depending on accomodation. Here we can write it
% out to a file to be read in later.
lensFile = sprintf('navarroAccomodated_%0.2f.dat',accommodation);
writeNavarroLensFile(navarroAccom,fullfile(workingFolder,lensFile));
fprintf('Wrote out a new lens file: \n')
fprintf('%s \n \n',fullfile(workingFolder,lensFile));

renderRecipe.camera.specfile.value = fullfile(workingFolder,lensFile);

end