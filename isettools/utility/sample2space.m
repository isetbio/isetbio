function [r, c] = sample2space(rSamples, cSamples, rowDelta, colDelta)
% Return the physical position of samples
%
% Syntax:
%   [r, c] = sample2space(rSamples, cSamples, rowDelta, colDelta)
%
% Description:
%    We treat the center of the samples as (0, 0) and use the sampling
%    spacing to calculate the locations of the other samples. The locations
%    r and c are in the same units as the deltas.
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

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting

% Examples:
%{
    oi = sceneCreate;
    [r, c] = sample2space(sceneGet(oi, 'rows'), sceneGet(oi, 'cols'), ...
      sceneGet(oi, 'hres'), sceneGet(oi, 'wres'));

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
