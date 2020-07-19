function [ior,wave] = legrandRefractiveIndices
% Calculate the IOR for the LeGrand model
%
% Synopsis
%   [ior,wave] = legrandRefractiveIndices
%
% Inputs
%
% Optional key/val pairs
%
% Outputs
%   ior:
%
% See also
%   legrandLensCreate; navarroLensCreate;

% Examples:
%{
[ior,wave] = legrandRefractiveIndices;
ieNewGraphWin; plot(wave,ior);
%}

%%
wave = (400:10:800); % nm

%% These are the index of refraction parameters

%           [Cornea; Aqueous; Lens; Vitreous]
n_inf    = [1.3610 1.3221 1.3999 1.3208]';
K        = [7.4147 7.0096 9.2492 6.9806]';
lambda_o = [130.0 130.0 130.0 130.0]';

% In the original code, but seems unused
%  V        = [56 53 50 53]';   

% Cornu dispersion equation
% Rows will be each ocular media
ior = n_inf + K ./ (wave - lambda_o);

% Put them in the columns, not rows
ior = ior';

%%