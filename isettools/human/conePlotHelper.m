function [xySubMosaic, coneSubMosaic, support,spread,delta] = conePlotHelper(sensor, conesToPlot, coneSize)
% Given a sensor struct, a desired number of cones to plot, and a desired
% cone size (diameter in pixels), return the variables expected as input
% arguments by the conePlot function
%
% Example: Generate an image of the central 30x30 cones, with each cone
% having a 20 pixel diameter
%
%   sensor = sensorCreate('human');
% 
%   conesToPlot = 30;
%   coneSize = 20;
% 
%   [xy,coneType, support,spread,delta] = conePlotHelper(sensor, conesToPlot, coneSize);
% 
%   whiteBackground = false;
%   [~,~,~, coneMosaicImage] = conePlot(xy,coneType, support,spread,delta, whiteBackground);
%   imshow(coneMosaicImage); truesize;

    sensorSeperation(1) = sensorGet(sensor, 'wspatial resolution', 'um');
    sensorSeperation(2) = sensorGet(sensor, 'hspatial resolution', 'um');
    xy = sensorGet(sensor, 'cone xy'); 
    xCoords = squeeze(xy(:,1));
    yCoords = squeeze(xy(:,2));
    selectConeIndices = find( (abs(xCoords) < sensorSeperation(1)*conesToPlot/2) & ...
                          (abs(yCoords) < sensorSeperation(2)*conesToPlot/2) ...
                            );   
    coneType = sensorGet(sensor, 'coneType');

    spread = 0.25*coneSize; support = round(spread*4*[1 1]); delta = 1.5/coneSize;
    coneSubMosaic = coneType(selectConeIndices);
    xySubMosaic = xy(selectConeIndices,:);
end