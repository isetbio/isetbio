function [macbethLAB, macbethXYZ, dE, vci] = ...
    macbethColorError(vci, illName, cornerPoints, method)
% Analyze color error of a MCC in the image processor window
%
% Syntax:
%	[macbethLAB, macbethXYZ, deltaE, vci] = ...
%        macbethColorError(vci, illName, cornerPoints, method)
%
% Description:
%    The user interactively identifies the position of the MCC in the image
%    processor window. This routine analyzes the RGB MCC data and plots the
%    results. The analyses are:
%         - Comparisons of the data and ideal CIELAB values
%         - A histogram of delta E values between the data and ideal
%         - Comparisons of the gray series L* values
%         - Comparisons of the chromaticity (xy) values
%
%    Examples are contained in the code. To access, type 'edit
%    macbethColorError.m' into the Command Window.
%
% Inputs:
%    vci          - (Optional) An image processor structure from ISET
%                   containing the processed MCC. ON the return, the corner
%                   points of the MCC are stored in the structure. Default
%                   is to retrieve a vcimage.
%    illName      - (Optional) The illuminant assumed for the MCC. Default
%                   is 'd65'.
%    cornerPoints - (Optional) A rect defining the outer points of the MCC.
%                   The order of the outer corners is white, gray, blue,
%                   brown. Default is to pull the mcc corner points from
%                   the vcimage.
%    method       - This defines the color space of the vci image data in
%                   the result field. Default 'sRGB'. You can set the
%                   method equal to 'custom' to use the monitor model
%                   stored in the vci.
%
% Outputs:
%    macbethLAB   - CIE LAB Errors
%    macbethXYZ   - XYZ Errors
%    dE           - delta E Errors
%    vci          - vci image errors
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - Could add display gamut to chromaticity plot?]
%    * TODO: FIX? vcimageMCCXYZ does not exist
%    * [Note: JNM - I have zero idea of what the outputs are supposed to
%      represent, so I am almost certain I've botched them.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    01/31/18  jnm  Formatting
%    04/07/18  dhb  ETTBSkip on broken example.

% Examples:
%{
    % ETTBSkip.  Example is broken. Function may or may not be broken too.
    vci = vcGetObject('vcimage');
	[macbethXYZ, whiteXYZ] = vcimageMCCXYZ(vci);
%}

%% Check variables
if notDefined('vci'), vci = vcGetObject('vcimage'); end
if notDefined('illName'), illName = 'd65'; end

% cornerPoints has the coordinates of the corners of the MCC.
if notDefined('cornerPoints')
    cornerPoints = imageGet(vci, 'mcc corner points');
end

% Either custom or sRGB calculation styles
if notDefined('method'), method = 'sRGB'; end

whiteIndex = 4;
gSeries = 4:4:24;

%% Retrieve ideal and image MCCdata LAB and XYZ values
% These are Processor window XYZ values using the monitor model. They are
% computed using a model monitor with peak luminance of Y = 1 cd/m2.
[macbethXYZ, whiteMacbethXYZ, cornerPoints] = vcimageMCCXYZ(vci, ...
    cornerPoints, method);
vci = imageSet(vci, 'mcc corner points', cornerPoints);

% The idealXYZ values of the MCC
idealXYZ = macbethIdealColor(illName, 'xyz');
%   macbethsRGB = xyz2srgb(reshape(idealXYZ, 4, 6, 3));
%   vcNewGraphWin; image(macbethsRGB);

% This is the white surface
whiteIdealXYZ = idealXYZ(whiteIndex, :);

% The max ideal luminance is whiteXYZ(2), so we scale by this quantity to
% put idealXYZ and macbethXYZ in the same range.
macbethXYZ = (macbethXYZ/whiteMacbethXYZ(2))*whiteIdealXYZ(2);
whiteMacbethXYZ = macbethXYZ(whiteIndex, :);
% vcNewGraphWin; plot(macbethXYZ(:), idealXYZ(:), '.')
% axis equal, grid on

% %% Initialize figure
% macbethEvaluationGraphs(L, rgb, idealRGB, sName);

vcNewGraphWin([], 'upperleft big'); clf;
set(gcf, 'name', sprintf('VCIMAGE: %s', imageGet(vci, 'name')))

%% LAB positions of the patches
subplot(2, 2, 1),
idealLAB = ieXYZ2LAB(idealXYZ, whiteIdealXYZ);
macbethLAB = ieXYZ2LAB(macbethXYZ, whiteMacbethXYZ);

plot(macbethLAB(:, 2), macbethLAB(:, 3), 'o')
line([idealLAB(:, 2), macbethLAB(:, 2)]', ...
    [idealLAB(:, 3), macbethLAB(:, 3)]');

xlabel('a (red-green)');
ylabel('b (blue-yellow)')
grid on;
axis square
title('CIELAB color plane')

%% Show histogram of delta E errors
subplot(2, 2, 2)

% We compute the delta E difference between the data and the ideal
dE = deltaEab(macbethXYZ, idealXYZ, whiteIdealXYZ);
hist(dE);
grid on;
axis square
title(sprintf('Mean deltaE = %.2f', mean(dE)));

%% Show the gray series L* values
subplot(2, 2, 3)
plot(1:6, macbethLAB(gSeries, 1), '-o', 1:6, idealLAB(gSeries, 1), 'x');
xlabel('Gray patch');
ylabel('L*');
axis square;
grid on;
title('Achromatic series')

%% Lines between chromaticities of the ideal and current Processor data
subplot(2, 2, 4)

% Exclude very black surfaces from chromaticity plot.
list = find(macbethXYZ(:, 2) > 0.01);
xy = chromaticity(macbethXYZ(list, :));
chromaticityPlot(xy, 'gray', 256, 0);

% Draw little lines to the ideal position
hold on;
idealxy = chromaticity(idealXYZ(list, :));
line([idealxy(:, 1), xy(:, 1)]', [idealxy(:, 2), xy(:, 2)]');

%% Store the data in the figure
uData.macbethXYZ =  macbethXYZ;
uData.macbethLAB = macbethLAB;
uData.idealLAB = idealLAB;
uData.deltaE = dE;

set(gcf, 'userdata', uData);

end