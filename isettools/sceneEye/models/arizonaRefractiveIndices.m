function [ior,wave,txt] = arizonaRefractiveIndices(accom)
% Return the IoR data for the Arizona model
%
% Synopsis
%
% Input
%   wave:   Wavelength sames, usually 400:10:800
%   accom:  Accomodation in diopters
%
% Optional key/val pairs
%
% Outputs
%    ior   - Matrix of length(wave), of the refractive indices
%            through the cornea, aqueous, lens, vitreous
%    wave   - Sample wavelengths (nm)
%    txt   - {'cornea','aqeous','lens','vitreous'};
%
% See also
%   navarroRefractiveIndices, legrandRefractiveIndices, arizonaWrite
%

% Examples"

%{
accom = 1
ior = arizonaRefractiveIndices(accom);
%}
%{
% Compare with other dispersion curves (plotting)

  legend('cornea', 'aqueous', 'lens', 'vitreous');
  xlabel('Wavelength (nm)'); ylabel('Index of Refraction');
  title('Dispersion Curves (Arizona vs Atchinson)')
  set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
  set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 2)
%}

%% The idea

% It is the same as in navarro and legrand cases

% We calculate the dispersion curves of each ocular media. In the lens
% file, each surface material is linked to an "ior slot" (ior1,
% ior2, etc.) When the ray is traveling through that material, it will
% follow the curve defined by the spectrum in the corresponding slot.

% Our convention (hard coded in writeArizonaLens) is always:
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous

% The Arizona eye changes the IOR of the lens with accom. The full
% dispersion curve is then calculated using the Abbe number. The full lens
% model can be found here:
% https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf

%%
wave = (400:10:800); % nm
txt  = {'cornea','aqeous','lens','vitreous'};

%%
% This is equivalent to IOR at 589.3 nm
% [cornea aqueous lens vitreous]
n_d = zeros(1, 4);
n_d(1) = 1.377;
n_d(2) = 1.337;
n_d(3) = 1.42 + (0.00256 * accom) - (0.00022 * accom ^ 2);
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
ior = zeros(numel(wave),4); 

for ii = 1:length(n_d)
    m = (n_d(ii)  - 1) / (V_d(ii) * (486.1 - 656.3));
    b = n_d(ii) - m * 589.3;
    ior(:,ii) = m * wave + b;

    % Compare with other dispersion curves (plotting)
    %{
    colors = {'r', 'b', 'm', 'g'};
    dispersionFig = ieNewGraphWin; 
    iorAtch = ieReadSpectra(...
    fullfile(piRootPath, 'data', 'lens', 'IORofEye.mat'), wave);
    
    plot(wave, ior, colors{ii}); 
    hold on; grid on;    
    plot(wave, iorAtch(:, ii), [colors{ii} '--'], ...
        'HandleVisibility', 'off');
    %}
end



end
