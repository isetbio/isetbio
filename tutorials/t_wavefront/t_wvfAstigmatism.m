function t_wvfAstigmatism
% Compute the wavefront-based PSF for various astigmatism and blur levels.
%
% Syntax:
%   t_wvfAstigmatism
%
% Description:
%    Compute the wavefront-based point-spread function for various
%    astigmatism and blur levels.
%
%    This illustrates the effect of Zernicke coefficients 4 and 5.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/12       (c) Wavefront Toolbox Team, 2012
%    01/01/18  dhb  Handled JNM notes.
%    09/27/18  jnm  Formatting

%% Initialize and set parameters
ieInit;

%% Range for plotting
maxUM  = 20;

%% Set up default parameters structure with diffraction limited default
% The ranges for coefficients here and below are reasonable given typical
% variation within human population.  If we look at the diagonal of the
% covariance matrix for coefficients that we get from the Thibos
% measurements (see wvfLoadThibosVirtualEyes we see that for the third
% through sixth coefficients, the standard deviations (sqrt of variances on
% the diagonal) range between about 0.25 and about 0.5.
wvfP = wvfCreate;
wvfParams = wvfComputePSF(wvfP);
z4 = -0.5:0.5:0.5;
z5 = -0.5:0.5:0.5;
[Z4, Z5] = meshgrid(z4, z5);
Zvals = [Z4(:), Z5(:)];

%% Alter defocus and astigmatism
% Make a plot of the psf for each case.
h = vcNewGraphWin;
set(h, 'Position', [0.5 0.5 0.45 0.45]);
wList = wvfGet(wvfParams, 'calc wave');
for ii = 1:size(Zvals, 1)
    wvfParams = wvfSet(wvfParams, 'zcoeffs', Zvals(ii, :), ...
        {'defocus' 'vertical_astigmatism'});
    wvfParams = wvfComputePSF(wvfParams);

    % Don't open a new window with each plot. Allow them to accumulate in
    % the subplots.
    subplot(3, 3, ii)
    wvfPlot(wvfParams, '2dpsfspace', 'um', wList, maxUM, 'nowindow');
    title(sprintf('Defocus = %.1f Astig == %.1f\n', Zvals(ii, 1), ...
        Zvals(ii, 2)));
end
