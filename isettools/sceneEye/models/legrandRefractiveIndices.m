function [ior,wave,txt] = legrandRefractiveIndices
% Calculate the IOR for the LeGrand model
%
% Synopsis
%   [ior,wave,txt] = legrandRefractiveIndices
%
% Inputs
%  N/A
%
% Optional key/val pairs
%  N/A
%
% Outputs
%   ior:   Index of refraction from 400:10:800
%   wave:  Wavelength samples
%   txt:   {'cornea','aqeous','lens','vitreous'}
%
% See also
%   legrandLensCreate; legrandWrte, navarroLensCreate;

% Examples:
%{
[ior,wave] = legrandRefractiveIndices;
ieNewGraphWin; plot(wave,ior);
%}

%%
wave = (400:10:800); % nm
txt  = {'cornea','aqeous','lens','vitreous'};

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

end