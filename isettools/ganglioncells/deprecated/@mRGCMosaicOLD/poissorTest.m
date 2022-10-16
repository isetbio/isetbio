
load('CronerKaplanNoiseMidget.mat');
CronerKaplanNoise = table2array(CronerKaplanNoise);

meanResponses = CronerKaplanNoise(:,1);
responseSTD = CronerKaplanNoise(:,2);
signalToNoiseRatioMidget =  meanResponses ./ responseSTD;

figure(1); clf;
hold on;
plot(1:100, sqrt(1:100), 'k--');
plot(meanResponses, responseSTD, 'ro'); 
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [1 100], 'YLim', [1 100], 'XTick', [1 3 10 30 100], 'YTick', [1 3 10 30 100], 'FontSize', 16);
axis 'square'
grid on


