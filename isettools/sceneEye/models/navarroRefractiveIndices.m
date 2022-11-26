function  [ior,wave,txt] = navarroRefractiveIndices(accom)
% Return the set of Navarro model index of refractions for an accommodation 
%
% Syntax:
%   [ior,wave,txt] = navarroRefractiveIndices(wave, accom)
%
% Description:
%    The the index of refraction in the Navarro model for the various
%    ocular media depends on the accomodative state.  This function
%    calculates the IoR for the four elements of the model.
%
% Inputs:
%    accom - Numeric. Non-navarro accommodation, in diopters.
%
% Optional key/value pairs:
%    None.
%
% Outputs:
%    ior   - Matrix of length(wave), of the refractive indices
%            through the cornea, aqueous, lens, vitreous
%    wave   - Sample wavelengths (nm)
%    txt   - {'cornea','aqeous','lens','vitreous'};
%
% See also
%   navarro*
%

% Examples:
%{
  % Is the dependence on accommodation is worth any effort?
  accDiopters = 1/10;
  [ior,wave] = navarroRefractiveIndices(accDiopters);
  ieNewGraphWin; plot(wave,ior); 

  accDiopters = 1/0.4;
  [ior,~,txt] = navarroRefractiveIndices(accDiopters);
  hold on; plot(wave,ior,'o'); 
  grid on;
  txt
%}

%% Fixed variables

wave = 400:10:800;
txt  = {'cornea','aqeous','lens','vitreous'};

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
if min(wave) < 350 || max(wave) > 800
    error(['Units for wavelength don''t seem correct...' ...
        ' they should be in nm roughly in the range 400 to 700.']);
end

% The refractive index equations use um, so let's convert to um. We want
% to stick with the nm convention though, to match ISET and PBRT.
waveUM = wave .* 10 ^ -3;

lambda0 = 0.1673;
L = 1 ./ (waveUM .^ 2 - lambda0 ^ 2);
RI_all = zeros(4, length(wave));
for i = 1:length(mediaNames)
    RI_all(i, :) = a(i, 1) + a(i, 2) * waveUM .^ 2 ...
        + a(i, 3) * L + a(i, 4) * L .^ 2;
end

% We used to return this.
%{
cor = RI_all(1, :);
aqu = RI_all(2, :);
len = RI_all(3, :);
vit = RI_all(4, :);
%}

ior = RI_all'; % [cornea(:),aqueuous(:), lens(:), vitreous(:)];

end
