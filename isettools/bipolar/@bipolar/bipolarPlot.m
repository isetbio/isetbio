function bipolarPlot(obj)
% Plot the input to and repsonse of a bipolar object.
% 
% 5/2016 JRG (c) isetbio team

vcNewGraphWin([],'upperLeftBig');

responseSize = size(obj.responseCenter);

responseRS = reshape(obj.responseCenter,responseSize(1)*responseSize(2),responseSize(3));
% responseRS = reshape(obj.responseSurround,responseSize(1)*responseSize(2),responseSize(3));
% responseRS = reshape(obj.responseCenter-obj.responseSurround,responseSize(1)*responseSize(2),responseSize(3));

plot(.001*(1:responseSize(3)),responseRS);
xlabel('Time (sec)'); 
ylabel('Response (AU)');
title('Bipolar Mosaic Response');