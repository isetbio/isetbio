function  [img, coords] = imageMakeMontage(hc, sliceList, nCols, backVal)
% Create an image comprising a montage of the slices in hc
%
% Syntax:
%   img = imageMakeMontage(hc, [sliceList], [nCols], [backVal])
%
% Description:
%    Create an image comprising a montage of the slices in hc
%
% Inputs:
%    hc        - x by y by sliceNum
%    sliceList - The number of slices to include. Default is all.
%    nCols     - number of columns in the image montage. Default is []
%    backVal   - Background color. Default is 0.
%
% Outputs:
%    img       - The image data
%    coords    - Image coordinate data
%
% Notes:
%    * [Note: JNM - Example isn't working due to no instantiated hypercube.
%      Can we get this fixed please? My macbeth chart addition allows
%      something to display, but I imagine that's not what is intended for
%      this example.]
%
% See also:
%   imageMontage
%

% History:
%    xx/xx/12       (c) Imageval, 2012
%    12/07/17  jnm  Formatting & edit example

% Examples:
%{
    % Example doesn't work without the hypercube being instantiated
    %
    % This is cheating since I'm using a Macbeth Chart instead of a
    % hypercube, but I only know how to instantiate the chart...
    hc = macbethChartCreate;
    [img, coords] = imageMakeMontage(hc.data);
	vcNewGraphWin;
    imagesc(img);
    colormap(gray);
    axis equal;
    axis off;
%}

%%
if notDefined('hc'), error('hypercube data required.'); end

[r, c, w] = size(hc);
if notDefined('sliceList'), sliceList = 1:w; end
if notDefined('nCols'), nCols = []; end
if notDefined('backVal'), backVal = 0; end

% If the hc is too large, refuse to play.
if(any(size(hc) > 10000))
    error('hc size too large - refusing to continue...');
end

%% Make a best guess about the number of columns
nImages = length(sliceList);
if(isempty(nCols)), nCols = ceil(sqrt(nImages) * sqrt((r / c))); end
nRows = ceil(nImages / nCols);

% Allocate for the image. Can we put the class(hc) in here?
img = ones(r*nRows, c*nCols, class(hc))*backVal;

count = 0;
nSlices = length(sliceList);
coords = zeros(nSlices, 2);
for ii = 1:nSlices
    curSlice = sliceList(ii);
    count = count + 1;
    x = rem(count - 1, nCols) * c;
    y = floor((count - 1) / nCols) * r;
    img(y + 1:y + r, x + 1:x + c) = hc(:, :, curSlice);
    coords(ii, :) = [x + 1, y + 1];
end

end