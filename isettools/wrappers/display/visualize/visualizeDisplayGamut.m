function visualizeDisplayGamut(primariesXYZ)
% Method to visualize the display's gamut
% Extract the maximum luminance for each primary (Y tristimulus value)
maxLuminanceCdPerM2 = primariesXYZ(:,2);
% Extract the (x,y) chromaticity coordinates, e.g. x = X / (X+Y+Z)
xChroma = primariesXYZ(:,1) ./ sum(primariesXYZ,2);
yChroma = primariesXYZ(:,2) ./ sum(primariesXYZ,2);

figure(); clf;
% Plot the (x,y) coords of the RGB guns
plot(xChroma(1), yChroma(1), 'rs', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]); hold on
plot(xChroma(2), yChroma(2), 'gs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 1 0.5]);
plot(xChroma(3), yChroma(3), 'bs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.5 1]);
% Plot the CIE (xy) background 
renderCIEdiagramBackground();
% Replot the (x,y) cords of the RGB guns
plot(xChroma(1), yChroma(1), 'rs', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0.5]); 
plot(xChroma(2), yChroma(2), 'gs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 1 0.5]);
plot(xChroma(3), yChroma(3), 'bs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.7 1]);
xx = [xChroma; xChroma(1)]; yy = [yChroma; yChroma(1)];
plot(xx,yy, 'k--', 'LineWidth', 1.5); hold on;
axis 'square'; axis 'xy';grid on; 
xlabel('\it x-chroma'); ylabel('\it y-chroma');
legend({...
    sprintf('R (max lum: %2.1f cd/m2)', maxLuminanceCdPerM2(1)) ...
    sprintf('G (max lum: %2.1f cd/m2)', maxLuminanceCdPerM2(2)) ...
    sprintf('B (max lum: %2.1f cd/m2)', maxLuminanceCdPerM2(3)) });
title('display gamut')
end

function renderCIEdiagramBackground()
% Method to render the shoe-horse CIE color background
wave = 420:5:700;
XYZcolorMatchingFunctions = ieReadSpectra('XYZ', wave);
xOutline = XYZcolorMatchingFunctions(:,1)./sum(XYZcolorMatchingFunctions,2);
yOutline = XYZcolorMatchingFunctions(:,2)./sum(XYZcolorMatchingFunctions,2);
xOutline(end+1) = xOutline(1);
yOutline(end+1) = yOutline(1);

N = 500;
x = (0:(N-1))/N;
[X,Y] = meshgrid(x);
[iCol, iRow] = meshgrid(1:N);
iCol = iCol(:); iRow = iRow(:);
X = X(:); Y = Y(:);
[in,on] = inpolygon(X(:),Y(:),xOutline,yOutline);
idx = find(in==true);
backgroundLuminance = 90;
lum = zeros(numel(idx),1) + backgroundLuminance/683;
XYZ = xyYToXYZ([X(idx) Y(idx) lum]');
c = xyz2rgb(XYZ');
c(c<0) = 0;
c(c>1) = 1;

backgroundImage = zeros(N,N,3)+1;
for ix = 1:numel(idx)
    if (in(idx(ix))) || (on(idx(ix)))
       theColor = c(ix,:);
       backgroundImage(iRow(idx(ix)),iCol(idx(ix)),:) = theColor;
    end
end
x = (0:(N-1))/N;
image(x,x,backgroundImage);
hold on
plot([0 0 0.9 0.9 0], [0 0.9 0.9 0 0], 'k-');
plot(xOutline,yOutline,'ko-', 'MarkerFaceColor', 'k');
indices = [1 10 14:2:35 37 57];
for k = 1:numel(indices)
    wI = indices(k);
    text(xOutline(wI)+0.01, yOutline(wI)+0.02, sprintf('%2.0f',wave(wI)));
end
set(gca, 'XLim', [0 0.9], 'YLim', [0 0.9], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1, 'FontSize', 16);
end
