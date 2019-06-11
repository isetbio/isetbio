function [cor, aqu, len, vit] = getNavarroRefractiveIndices(wave, accom)
% Return the Navarro IoR at given accommodation and wavelengths
%
% Syntax:
%   [cor, aqu, len, vit] = getNavarroRefractiveIndices(wave, accom)
%
% Description:
%    Returns the index of refraction (from Navarro's model) for the various
%    ocular media at the given wavelengths and the given accomodative state
%
% Inputs:
%    wave  - Array. A wavelengths vector in nm.
%    accom - Numeric. Non-navarro accommodation, in diopters.
%
% Outputs:
%    cor   - Array. An array of length(wave), of the refractive indices
%            through the cornea.
%    aqu   - Array. An array of length(wave), of the refractive indices
%            through the aqueous solution.
%    len   - Array. An array of length(wave), of the refractive indices
%            through the lens.
%    vit   - Array. An array of length(wave), of the refractive indices
%            through the vitreous fluid.
%
% Optional key/value pairs:
%    None.
%

%% Convert from Herzberger
% Let's start by converting from the Herzberger coefficients described in
% Navarro et al. 1985 into the coefficients used by Zemax.

% From Atchinson 2005
% a_n = A0 + A1 * lambda ^ 2 + P * L + R * L ^ 2
% [A0 A1 P R] (col) for a_1, a_2, a_3, a_4 (row)
A = [ 0.66147196 -0.40352796 -0.28046790  0.03385979
     -4.20146383  2.73508956  1.50543784 -0.11593235
      6.29834237 -4.69409935 -1.57508650  0.10293038
     -1.75835059  2.36253794  0.35011657 -0.02085782];

n = [ 1.3975 1.3807  1.37405 1.3668
      1.3593 1.3422  1.3354  1.3278
      1.4492 1.42625 1.4175  1.4097
      1.3565 1.3407  1.3341  1.3273];

% Add effect of accommodation, if necessary
n(3, :) = n(3, :) + 9e-5 * (10 * accom + accom ^ 2);

% For each ocular media
mediaNames = {'cornea', 'aqueous', 'lens', 'vitreous'};
% Rows are ocular media, columns are coefficients (e.g. A, D, B, C)
a = zeros(4, 4);
for i = 1:4
    % [n**, n_F, n_C, n*]
    n_curr = n(i, :);

    % [A D B C]
    a(i, :) = [n_curr * (A(:, 1)), n_curr * (A(:, 2)), ...
        n_curr * (A(:, 3)), n_curr * (A(:, 4))];
%      fprintf('%s:    A = %f    B = %f    C = %f    D = %f    \n', ...
%          mediaNames{i}, a(i, 1), a(i, 3), a(i, 4), a(i, 2));
end

%% Calculate refractive indices at given wavelength
% Check units for wavelength
if(sum(wave < 400 + wave > 800) > 0)
    error(['Units for wavelength don''t seem correct...' ...
        ' they should be nm here.']);
end

% The equations use um, so let's convert to um. We want to stick with the
% nm convention though, to match ISET and PBRT.
wave = wave .* 10 ^ -3;

lambda0 = 0.1673;
L = 1 ./ (wave .^ 2 - lambda0 ^ 2);
RI_all = zeros(4, length(wave));
for i = 1:length(mediaNames)
    RI_all(i, :) = a(i, 1) + a(i, 2) * wave .^ 2 ...
        + a(i, 3) * L + a(i, 4) * L .^ 2;
end

cor = RI_all(1, :);
aqu = RI_all(2, :);
len = RI_all(3, :);
vit = RI_all(4, :);

end
