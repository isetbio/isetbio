function [figH, m, cbH] = imageMontage(hc, slices, numCols, figNum)
% Create a window with a montage of the slices in the hypercube data
%
% Syntax:
%   [figH, m, cbH] = imageMontage(hc, [slices], [numCols], [figNum])
%
% Description:
%    Create a window with a montage of the slices in the hypercube data
%
%    Examples are located within the code. To access the examples, type
%    'edit imageMontage.m' into the Command Window.
%
% Inputs:
%    hc      - Hypercube data
%    slices  - (Optional) Slices to create a montage with. Default is [].
%    numCols - (Optional) Number of columns in the montage. Default is [].
%    figNum  - (Optional) Specify the figure. Default is Figure's Default.
%
% Outputs:
%    figH    - The hypercube figure
%    m       - The image montage
%    cbH     - The colorbar for the hypercube
%
% Optional key/value pairs:
%    None
%
% Notes:
%    * [Note: JNM - removed listed input variables that were not actually
%      in the header. Variables removed: wavebands, cmap, crop, flip.]
%
% See Also:
%    imageMakeMontage
%

% History:
%    xx/xx/12       (c) Imageval, 2012
%    12/07/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match the Wiki.

% Examples:
%{
    hc = macbethChartCreate;
    imageMontage(hc.data);
    colormap(gray);
    axis image
%}

if notDefined('slices'), slices = []; end
if ~exist('numCols', 'var'), numCols = []; end
if ~exist('figNum', 'var'), figH = figure; else, figH = figure(figNum); end

%% 
m = imageMakeMontage(hc, slices, [], numCols);
imagesc(double(m));
axis image;
cbH = colorbar;

end