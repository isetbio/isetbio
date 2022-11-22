function thisR = setNavarroAccommodation(thisR, accommodation, workingFolder)
% Change renderRecipe to match the accommodation
%
% Syntax:
%   thisR = setNavarroAccommodation(thisR, accommodation, [workingFolder])
%
% Description:
%    For the human eye models, when accommodation changes the lens
%    file changes. We write these new files out and reference them in
%    the structure. 
%
%   This scope of the Navarro model includes accommodation from 0 to
%   10 diopters. A value of 0 diopters means the focus is at infinity.
%   10 diopters means the eye model focus is at 0.1 meter. 
%
%   However, the Navarro accommodation values do not match the Zemax
%   values.  So we call a function that converts the user's
%   accommodation value to the one from Zemax.  We should probably
%   allow ourselves to turn this conversion off (BW).
%
% Inputs:
%    thisR         - Object (Render recipe)
%    accommodation - Numeric. The accommodation to shape the modified
%                    thisR by.
%
% Optional
%    workingFolder - String. By default it is
%                    thisR.get('lens output dir')
%
% Outputs:
%    thisR          - The modified render recipe.
%
%
% See also
%   navarroLensCreate

% History:
%    xx/xx/17  TL   Created by Trisha Lian IESTBIO Team 2017
%    12/19/17  jnm  Formatting
%    01/25/19  JNM  Changed output of filenames to remove extra periods,
%                   which is breaking Windows executions.
%    05/29/19  JNM  Second documentation pass (minor tweaks)

%% Check and make sure this recipe has a human eye model
if ~strcmp(thisR.get('camera subtype'),'humaneye')
    warning('The camera type is not a human eye model. Returning untouched.');
    return;
end

%% Check inputs
if(~(accommodation >= 0 && accommodation <= 10))
    % This is the scope of the model.  0 diopters means the focus is at
    % infinity.  10 diopters means the eye model focus is at 0.1 meter.
    error('Accommodation must be between 0 and 10 diopters.');
end

% Default output directory
if notDefined('workingFolder'), workingFolder = thisR.get('lens output dir'); end
if ~exist(workingFolder, 'dir')
    error('Working folder does not exist.');
end

% Maybe set a switch to include this or not.
accommodation = convertToNavarroAccomm(accommodation);

% Simply over-write the eye model with the new accommodation
navarroWrite(thisR,accommodation);

end

% The original code that has now been replaced.
%
%{
%% Convert accommodation

% See the function description for more information on why this is needed.
% Basically, the Navarro units are not a simple match to the focal distance
% when we compared them to the values computed by Zemax.
accommodation = convertToNavarroAccomm(accommodation);

%% Write out ocular media spectra files

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
wave = (400:10:800);    % nm
[cor, aqu, len, vit] = getNavarroRefractiveIndices(wave, accommodation);

% accommodation is 1/focal distance and thus has units of 1/m.
% 
% We represent the accommodation into the file name in this format
%
%     ior1_X_YY.spd 
%
% where X_YY becomes a value of X.YY and a focal distance of 1/(X.YY)
accStr = sprintf('%0.2f', accommodation);
accStr = strrep(accStr, '.', '_');
iorNames = {sprintf('ior1_%sdp.spd', accStr), ...
    sprintf('ior2_%sdp.spd', accStr), ...
    sprintf('ior3_%sdp.spd', accStr), ...
    sprintf('ior4_%sdp.spd', accStr)};

% I think rtb is an older reference to Render Toolbox
rtbWriteSpectrumFile(wave, cor, fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, aqu, fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, len, fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, vit, fullfile(workingFolder, iorNames{4}));

thisR.camera.ior1.value = fullfile(workingFolder, iorNames{1});
thisR.camera.ior2.value = fullfile(workingFolder, iorNames{2});
thisR.camera.ior3.value = fullfile(workingFolder, iorNames{3});
thisR.camera.ior4.value = fullfile(workingFolder, iorNames{4});

%% Attach lens file and set retina radius

% For Navarro (and Arizona), the lens file will change depending on
% accomodation. Here we can write it out to a file to be read in later.
% lensFile = sprintf('navarroAccomodated_%s.dat', accStr);
lensFile = navarroWrite(thisR,accommodation);

% writeNavarroLensFile(navarroAccom, fullfile(workingFolder, lensFile));
fprintf('Wrote out a new lens file: %s: \n',lensFile)
% fprintf('%s \n \n', fullfile(workingFolder, lensFile));

thisR.camera.lensfile.value = lensFile;
thisR.camera.lensfile.type = 'string';

%}

