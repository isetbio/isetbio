function visualizeDisplayLUT(gammaTable, direction)
% Method to visualize the display's LUT
figure(); clf;
% generate linearRGBvalues
lutLength = size(gammaTable,1);
linearRGBvalue = (0:(lutLength-1))/lutLength;
if (strcmp(direction, 'inverse'))
    gammaTable = gammaTable/max(gammaTable(:));
end
% plot the RGB LUTs
plot(linearRGBvalue, gammaTable(:,1), 'r-', 'LineWidth', 1.5); hold on;
plot(linearRGBvalue, gammaTable(:,2), 'g-', 'LineWidth', 1.5);
plot(linearRGBvalue, gammaTable(:,3), 'b-', 'LineWidth', 1.5);
grid on
% label plot
if (strcmp(direction, 'forward'))
xlabel('\it RGB (settings)'); 
ylabel('\it rgb (primaries)'); 
else
ylabel('\it RGB (settings)'); 
xlabel('\it rgb (primaries)'); 
end
axis 'square'
set(gca, 'FontSize', 16);
legend({'R primary', 'G primary', 'B primary'}, 'Location', 'NorthWest');
title(sprintf('%s gamma', direction))
end