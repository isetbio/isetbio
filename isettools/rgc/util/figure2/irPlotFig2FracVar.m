function irPlotFig2FracVar(experimentI,cellTypeI,fractionalVariance)

vcNewGraphWin([],'upperleftbig'); 
hold on;

badInds = fractionalVariance{experimentI,2,cellTypeI}<0;
fractionalVariance{experimentI,2,cellTypeI}(badInds) = 0;

switch cellTypeI
    case 1
        scatter(fractionalVariance{experimentI,1,cellTypeI},fractionalVariance{experimentI,2,cellTypeI},'r');
        scatter(fractionalVariance{experimentI,1,cellTypeI},fractionalVariance{experimentI,2,cellTypeI},8,'b');
    case 2
        scatter(fractionalVariance{experimentI,1,cellTypeI},fractionalVariance{experimentI,2,cellTypeI},'g','filled');
        scatter(fractionalVariance{experimentI,1,cellTypeI},fractionalVariance{experimentI,2,cellTypeI},8,'k','filled');
end

plot(0:.1:1,0:.1:1);
axis([0 1 0 1]);
xlabel('WN Score (AU)'); ylabel('NSEM Score (AU)');
title('Fractional Variance');
set(gca,'fontsize',16);