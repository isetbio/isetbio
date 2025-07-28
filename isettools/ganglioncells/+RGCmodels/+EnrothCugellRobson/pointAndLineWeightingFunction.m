function [pWFcenter, pWFsurround, lWFcenter, lWFsurround, STF, rDegs] = pointAndLineWeightingFunction(xDegs, yDegs, kc, ks, rcDegs, rsDegs, sfSPD, generatePlot)

    if (isempty(xDegs) && isempty(yDegs))
        error('either xDegs or yDegs must be non-empty');
    elseif (~isempty(xDegs) && isempty(yDegs))
        yDegs = xDegs;
    elseif (isempty(xDegs) && ~isempty(yDegs))
        xDegs = yDegs;
    end

    [x,y] = meshgrid(xDegs, yDegs);
    rDegs = sqrt(x.^2 + y.^2);
    pWFcenter = kc * exp(-(rDegs/rcDegs).^2);
    pWFsurround = ks * exp(-(rDegs/rsDegs).^2);

    lWFcenter = kc * rcDegs * sqrt(pi) * exp(-(xDegs/rcDegs).^2);
    lWFsurround = ks * rsDegs * sqrt(pi) * exp(-(xDegs/rsDegs).^2);

    STF = kc * rcDegs^2 * exp( - (pi * rcDegs * sfSPD).^2 ) - ks * rsDegs^2 * exp( - (pi * rsDegs * sfSPD).^2 );

    if (generatePlot)
        [~,idx] = max(pWF(:));
        [row,col] = ind2sub(size(pWF), idx);

        figure(1); clf;
        subplot(1,3,1);

        imagesc(xDegs, yDegs, pWF);
        hold on;
        plot(xDegs, yDegs(row) + xDegs*0, 'k-');
        plot(xDegs(col)+yDegs*0, yDegs, 'k-');
        axis 'image';
        title('point weighting function');

        subplot(1,3,2);

        plot(xDegs, pWF(row,:), 'k-', 'LineWidth', 1.5);
        title('point weighting function slice')


        subplot(1,3,3);
        dy = yDegs(2)-yDegs(1);
        plot(xDegs, lWFcenter, 'r-', 'LineWidth', 2); hold on;
        plot(xDegs, -lWFsurround, 'b-', 'LineWidth', 2);
        plot(xDegs, lWFcenter-lWFsurround, 'k-', 'LineWidth', 3.0);
        plot(xDegs, sum(pWF,1)*dy, 'm:', 'LineWidth', 1.5);
        legend({'lWFcenter', 'lWFsurround', 'lWFcenter - lWFsurround', 'sum(pWF,1)'})
        title('line weighting functions');
    end




end
