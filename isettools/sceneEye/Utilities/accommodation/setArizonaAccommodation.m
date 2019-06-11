function renderRecipe = setArizonaAccommodation(...
    renderRecipe, accommodation, workingFolder)
% Change renderRecipe to match the accommodation for the Arizona eye
%
% Syntax:
%   renderRecipe = setArizonaAccommodation(..
%       renderRecipe, accommodation, workingFolder)
%
% Description:
%    We change the fields of the renderRecipe to match accommodation. As
%    accommodation changes, the lens file will change, as will the index of
%    refraction for the lens media. We write these new files out and
%    reference them in the structure.
%
% Inputs:
%    renderRecipe  - Object. The un-modified renderRecipe.
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
%    05/29/19  JNM  Documentation pass

%% Check and make sure this recipe includes a realisticEye
if(~strcmp(renderRecipe.camera.subtype, 'realisticEye'))
    warning('The camera type is not a realisticEye. Returning untouched.');
    return;
end

%% Check inputs
if ~exist(workingFolder, 'dir')
    error('Working folder does not exist.');
end

%% Write out ocular media spectra files
% We calculate the dispersion curves of each ocular media. In the lens
% file, each surface material is linked to an "ior slot" (ior1,
% ior2, etc.) When the ray is traveling through that material, it will
% follow the curve defined by the spectrum in the corresponding slot.

% Our convention (hard coded in writeArizonaLens) is always:
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous

% The Arizona eye changes the IOR of the lens with accommodation. The full
% dispersion curve is then calculated using the Abbe number. The full lens
% model can be found here:
% https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf

wave = (400:10:800); % nm

% This is equivalent to IOR at 589.3 nm
% [cornea aqueous lens vitreous]
n_d = zeros(1, 4);
n_d(1) = 1.377;
n_d(2) = 1.337;
n_d(3) = 1.42 + (0.00256 * accommodation) - (0.00022 * accommodation ^ 2);
n_d(4) = 1.336;

% Abbe number
% [cornea aqueous lens vitreous]
V_d = zeros(1, 4);
V_d(1) = 57.1;
V_d(2) = 61.3;
V_d(3) = 51.9;
V_d(4) = 61.1;

% Calculate dispersion curves
% n_f is IOR at 486.1 nm
% n_c is IOR at 656.3 nm
% V_d  = (n_d - 1)/(n_f - n_c)
% V_d = ( f(589.3) - 1 )/( f(486.1) - f(656.3) )
ior = cell(1, 4); 

% Compare with other dispersion curves (plotting)
%{
colors = {'r', 'b', 'm', 'g'};
dispersionFig = figure(); hold on; grid on;
iorAtch = ieReadSpectra(...
    fullfile(piRootPath, 'data', 'lens', 'IORofEye.mat'), wave);
%}

for ii = 1:length(n_d)
    m = (n_d(ii)  - 1) / (V_d(ii) * (486.1 - 656.3));
    b = n_d(ii) - m * 589.3;
    ior{ii} = m * wave + b;

    % Compare with other dispersion curves (plotting)
    %{
    figure(dispersionFig);
    plot(wave, ior{ii}, colors{ii});  
    plot(wave, iorAtch(:, ii), [colors{ii} '--'], ...
        'HandleVisibility', 'off');
    %}
end

% Compare with other dispersion curves (plotting)
%{
legend('cornea', 'aqueous', 'lens', 'vitreous');
xlabel('Wavelength (nm)'); ylabel('Index of Refraction');
title('Dispersion Curves (Arizona vs Atchinson)')

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)
%}

iorNames = {sprintf('ior1_%0.2fdp_arizona.spd', accommodation), ...
    sprintf('ior2_%0.2fdp_arizona.spd', accommodation), ...
    sprintf('ior3_%0.2fdp_arizona.spd', accommodation), ...
    sprintf('ior4_%0.2fdp_arizona.spd', accommodation)};

rtbWriteSpectrumFile(wave, ior{1}, fullfile(workingFolder, iorNames{1}));
rtbWriteSpectrumFile(wave, ior{2}, fullfile(workingFolder, iorNames{2}));
rtbWriteSpectrumFile(wave, ior{3}, fullfile(workingFolder, iorNames{3}));
rtbWriteSpectrumFile(wave, ior{4}, fullfile(workingFolder, iorNames{4}));

renderRecipe.camera.ior1.value = fullfile(workingFolder, iorNames{1});
renderRecipe.camera.ior2.value = fullfile(workingFolder, iorNames{2});
renderRecipe.camera.ior3.value = fullfile(workingFolder, iorNames{3});
renderRecipe.camera.ior4.value = fullfile(workingFolder, iorNames{4});

%% Attach lens file 
% The lens file will change depending on accomodation. Here we can write it
% out to a file to be read in later.
lensFile = sprintf('arizonaAccomodated_%0.2f.dat', accommodation);
writeArizonaLensFile(accommodation, fullfile(workingFolder, lensFile));
fprintf('Wrote out a new lens file: \n')
fprintf('%s \n \n', fullfile(workingFolder, lensFile));

renderRecipe.camera.lensfile.value = fullfile(workingFolder, lensFile);
renderRecipe.camera.lensfile.type = 'string';

end