function [yNoise, mRGB] = macbethLuminanceNoise(vci, pointLoc)
% Analyze luminance noise in gray series of MCC from image processor window
%
% Syntax:
%	[yNoise, mRGB] = macbethLuminanceNoise([vci], [pointLoc])
%
% Description:
%    Analyze the luminance noise in gray series of MCC from the image
%    processor window
%
% Inputs:
%    vci      - (Optional) Image process structure. Default is to create a
%               new one using vcGetObject.
%    pointLoc - (Optional) Macbeth point locations. Default is [].
% 
% Outputs:
%	 yNoise   - Luminance noise
%    mRGB     - Linear RGB values of the display
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - Could add display gamut to chromaticity plot?]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    01/30/18  jnm  Formatting


%% Arguments
if notDefined('vci'), vci = vcGetObject('vcimage'); end
if notDefined('pointLoc'), pointLoc = []; end

%% Return the full data from all the patches
mRGB = macbethSelect(vci, 0, 1, pointLoc);

% Compute the standard deviation and mean for each patch. The ratio is the
% contrast noise.
gSeries = 4:4:24;
yNoise = zeros(1, length(gSeries));
for ii = 1 : length(gSeries)
    rgb = mRGB{gSeries(ii)};
    % Convert linear RGB values of the display to XYZ and then luminance
    macbethXYZ = imageRGB2XYZ(vci, rgb);
    Y = macbethXYZ(:, 2);
    
    % Calculate noise
    yNoise(ii) = 100 * (std(Y) / mean(Y));
end

%% Show it
vcNewGraphWin;
str = sprintf('%s: MCC luminance noise', imageGet(vci, 'name'));
set(gcf, 'name', str);
plot(yNoise);
line([1 6], [3 3], 'Linestyle', '--')
grid on
xlabel('Gray patch (white to black)')
ylabel('Percent luminance noise (std(Y)/mean(Y))x100');
legend({'data', '1000 photon (33 db) '}, 'Location', 'NorthWest')

end