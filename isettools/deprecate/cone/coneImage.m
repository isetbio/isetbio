function [uData, hf] = coneImage(obj, hf)
error('This function require some more debugging, do not use it');

% we first plot to a hidden figure
hiddenF = figure('visible', 'off');

% get position of the cones in um
x = obj.coneLocs(:, 1) * 1e6; y = obj.coneLocs(:, 2) * 1e6;
coneWidth = obj.pigment.width * 1e6;
coneHeight = obj.pigment.height * 1e6;
pdWidth = obj.pigment.pdWidth * 1e6;

% set image to high resolution
curUnits = get(gca, 'Units'); set(gca, 'Units', 'Points');
set(hiddenF, 'Position', [0 0 1200 1200]);
pos = [100 100 1000 1000];

% adjust aspect ratio
if pos(4) < pos(3) * obj.height / obj.width;
    pos(3) = pos(4) * obj.width / obj.height;
else
    pos(4) = pos(3) * obj.height / obj.width;
end

% compute conversion factor between meters and points
set(gca, 'Position', pos); set(gca, 'Units', curUnits);
meter2point = @(xm) xm * pos(3) / (max(x) - min(x) + coneWidth);

% scatter plot the cones
coneColor = [0 0 0; 0.9 0.1 0.1; 0.1 0.8 0.2; 0.1 0.1 1];
uData.s = scatter(x, y, meter2point(pdWidth)^2, ...
    coneColor(obj.pattern, :), 'fill', 'o');

% set limits of x/y axis based on the position of cones
xlim([min(x) - coneWidth/2,  max(x) + coneWidth/2]);
ylim([min(y) - coneHeight/2, max(y) + coneHeight/2]);
axis equal; xlabel('Position (um)'); ylabel('Position (um)');
set(gca, 'Color', 'black');
set(gca, 'FontSize', 20);

% get image
set(gca, 'Units', 'Pixels'); ti = get(gca, 'TightInset');
rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
set(gca, 'Units', curUnits);
frame = getframe(gca, rect);
uData.mosaicImage = frame.cdata;

% show image if hf is given
if isgraphics(hf, 'figure'), figure(hf); imshow(frame.cdata);
elseif isgraphics(hf, 'axes'), axes(hf); imshow(frame.cdata);
end
end