function [r, c] = sample2space(rSamples, cSamples, rowDelta, colDelta)
% Return the physical position of samples
%
% Syntax:
%   [r, c] = sample2space(rSamples, cSamples, rowDelta, colDelta)
%
% Description:
%    We treat the center of the samples as (0, 0) and use the sample
%    spacing to calculate the locations of the other samples. The locations
%    r and c are in the same units as the deltas.
%
%    For rows and columns separately, this routine uses the mean of the
%    passed sample values as the center location, and subtracts this from
%    the samples.  It then multiples by the delta. 
%
%    The use of the mean is fine, but note that it differs a little from
%    Matlab's conventions for the center of an array, which wants
%    coordinate 0 at location floor(N/2)+1, where N is the dimension of the
%    spatial support. At least, this is Matlab's convention for the FFT,
%    prior to application of ifftshift.
%
%    Examples in code.
%
% Inputs:
%    rSamples - Row Samples
%    cSamples - Column Samples
%    rowDelta - Row Delta
%    colDelta - Column Delta
%
% Outputs:
%    r        - Row spatial data
%    c        - Column spatial data
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting
%    01/16/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    [r, c] = sample2space(...
      1:sceneGet(scene, 'rows'), 1:sceneGet(scene, 'cols'), ...
      sceneGet(scene, 'hres'), sceneGet(scene, 'wres'));

    % Note unfortunate terminology in sceneGet:
    %   hres is height resolution; 
    %   wres is width resolution

    % Calculate the positions of every pixel this way:
    [X, Y] = meshgrid(c, r);
%}
%{
    cSamples = [1:.2:64];
    rSamples = [1:.2:64];
    rowDeltaMicrons = 5;
    colDeltaMicrons = 5;
    [rMicrons, cMicrons] = sample2space(rSamples, cSamples, ...
        rowDeltaMicrons, colDeltaMicrons)
%}

% Find the center.
%
% Note that although this method 
%   rCenter = max(rSamples(:)) / 2; 
%   cCenter = max(cSamples(:)) / 2;
% is faster than what we do, our method is correct even if the positions
% start at 1
rCenter = mean(rSamples(:));
cCenter = mean(cSamples(:));

% Convert to microns relative to center.
r = (rSamples - rCenter) * rowDelta;
c = (cSamples - cCenter) * colDelta;

end
