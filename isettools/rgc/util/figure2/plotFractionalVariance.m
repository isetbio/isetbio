function plotFractionalVariance(experimentI,cellTypeI,fractionalVariance)

vcNewGraphWin([],'upperleftbig'); 
hold on;

switch cellTypeI
    case 1
        scatter(fractionalVariance{experimentI,cellTypeI,1},fractionalVariance{experimentI,cellTypeI,2},'r');
        scatter(fractionalVariance{experimentI,cellTypeI,1},fractionalVariance{experimentI,cellTypeI,2},8,'b');
    case 2
        scatter(fractionalVariance{experimentI,cellTypeI,1},fractionalVariance{experimentI,cellTypeI,2},'g','filled');
        scatter(fractionalVariance{experimentI,cellTypeI,1},fractionalVariance{experimentI,cellTypeI,2},8,'k','filled');
end

plot(0:.1:1,0:.1:1);
axis([0 1 0 1]);
xlabel('WN Score (AU)'); ylabel('NSEM Score (AU)');
title('Fractional Variance');
set(gca,'fontsize',16);