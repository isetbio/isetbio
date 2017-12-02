function [row, col] = space2sample(rMicrons, cMicrons, pixelHeight, ...
    pixelWidth)
% (OBSOLETE) Convert spatial position (microns) into sample position
%
% Syntax:
%   [row, col] = space2sample(rMicrons, cMicrons, pixelHeight, pixelWidth)
%
% Description:
%    (OBSOLETE) This routine inverts sample2space. It converts spatial
%    positioning (in microns) into a sample position.
%
% Inputs:
%    rMicrons    - Row(s) of positions, in microns
%    cMicrons    - Column(s) of positions, in microns
%    pixelHeight - Height of pixels
%    pixelWidth  - Wifth of pixels
%
% Outputs:
%    row         - Row(s) of sample positions
%    col         - Column(s) of sample positions
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting

% Examples:
%{
    oi = sceneCreate;
    [row, col] = ...
        space2sample([1:10], [1:10], ...
            sceneGet(oi, 'hres'), sceneGet(oi, 'wres'));
	[X, Y] = meshgrid(row, col)
%}

disp('obsolete');

% tmp = (1 - 1 / 2 + rMicrons / pixelHeight); 
% row = tmp + max(tmp(:));
% tmp = (1 - 1 / 2 + cMicrons / pixelWidth);  
% col = tmp + max(tmp(:));

tmp = rMicrons / pixelHeight; % rescaling
row = tmp - tmp(1); % assuming samples are starting at 0;
tmp = cMicrons / pixelWidth;
col = tmp - tmp(1);

return;
