function visualizeWithEyeMovementSequence(obj, emTimeAxis)
% Visualize the eye movements
%
% Syntax:
%   visualizeWithEyeMovementSequence(obj, emTimeAxis)
%
% Description:
%    Visualize using the eye movement sequence
%
% Inputs:
%    obj        - Object. The oi Sequence object.
%    emTimeAxis - Vector. The eye movement time axis.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

hFig = figure();
set(hFig, 'Color', [1 1 1], 'Position', [10 10 1500 200]);
subplot('Position', [0.01 0.10 0.98 0.8]);

dtOI = 1000 * (obj.timeAxis(2) - obj.timeAxis(1));
maxX = 0;

hold on;
for oiIndex = 1:obj.length
    x = 1000 * obj.timeAxis(oiIndex) + [0 dtOI];
    if (x(2) > maxX), maxX = x(2); end
    y = [-0.5 0.5];
    plot([x(1) x(1) x(2) x(2) x(1)], [y(1) y(2) y(2) y(1) y(1)], 'k-', ...
        'LineWidth', 1.5);
end

dtEM = 1000 * (emTimeAxis(2) - emTimeAxis(1));

for emIndex = 1:numel(emTimeAxis)
    x = 1000 * emTimeAxis(emIndex) + [0 dtEM];
    if (x(2) > maxX), maxX = x(2); end
    y = [-1 1];
    plot([x(1) x(1) x(2) x(2) x(1)], [y(1) y(2) y(2) y(1) y(1)], 'r-', ...
        'LineWidth', 1.5);
end

dtMin = min([dtOI dtEM]);
box on;
grid on
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XLim', 1000 * [min([min(emTimeAxis) min(obj.timeAxis)]) - ...
    1 / 1000 maxX / 1000 + 1 / 1000], 'YLim', [-1.5 1.5]);
set(gca, 'FontSize', 14);
title(sprintf(['\\fontsize{16} \\color{red}eye movement sequence ' ...
    '(n = %d, integrationTime: %2.1f ms) \\color{black}; oi ' ...
    'sequence (n = %d, stim. interval: %2.1f ms)'], ...
    numel(emTimeAxis), dtEM, numel(obj.timeAxis), dtOI))

end
