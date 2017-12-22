function [res, cmap] = imagescOPP(oppImg, gam, nTable)
% NOT FINISHED IMPLEMENTING. Display three opponent-colors images 
%
% Syntax:
%   [res, cmap] = imagescOPP(oppImg, [gam], [nTable])
%
% Description:
%    Not yet fully implemented ...
%    
%    Create the three color-maps (black and white, blue and yellow, red and
%    green) for the provided image.
%
% Inputs:
%    oppImg   - The image to display in the opponent-colors images
%    gam      - (Optional) The Gamma exponent (Luminance). Default is 0.3.
%    nTable   - (Optional) The number of elements in the color map/table.
%               Default is 256.
%
% Outputs:
%    res      - The image
%    cmap     - The image color map information
%
% Notes:
%    * [Note: JNM - The example contained an example for imagescRGB, but
%      there was no reference to this function itself. I have since
%      rectified that.]
%    * [Note: JNM - The text in the description was for imagescRGB, as was
%      the example. I hope I have caught and removed all of the
%      inconsistencies and language related to the other function.]
%    * [Note: JNM - Is the 'not fully implemented ...' notes resolved? I've
%      left them for the time being, but I would like to remove them.]
%
% See Also:
%    imagescRGB

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010.
%    12/08/17  jnm  Formatting & rewrite description text

% Examples:
%{
    foo = load('trees');
    [r, c] = size(foo.X);
    for ii = 1:3, rgb(:, :, ii) = reshape(foo.map(foo.X, ii), r, c); end

    [imOp opMap] = imagescOPP(rgb);
    [imOp opMap] = imagescOPP(rgb, 1);
%}


if notDefined('oppImg'), error('Opponent image required.'); end
if notDefined('gam'), gam = 0.3; end
if notDefined('nTable'), nTable = 256; end

res = zeros(size(oppImg));
cmap = zeros(nTable, 3, 3);

% Create the images and the maps Three color maps (cm) are built. One for
% black-white (luminance), and others for red/green and blue/yellow. 
cm = {'bw', 'rg', 'by'};

for ii = 1:3
    % Color map with a gamma specified.
    cmap(:, :, ii) = ieCmap(cm{ii}, nTable, gam);

    tmp = oppImg(:, :, ii);
    mx  = max(abs(tmp(:)));
    
    % The luminance image doesn't get scaled much.
    if ii == 1
        res(:, :, ii)  = (tmp / mx) * nTable;
    else
        % Scale the opponent images so that 0 maps to 0.5, the range is
        % between -1 and 1.
        res(:, :, ii)  = (0.5 * (tmp / mx) + 0.5) * nTable;
    end
    
    vcNewGraphWin; 
    image(res(:, :, ii)); 
    colormap(cmap(:, :, ii));
    axis image;
    axis off;
    truesize
end

end