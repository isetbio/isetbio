function [hFig4, hFig5] = plotRFinstances(obj, eccRange)

   % Retrieve the data
   eccDegs = obj.synthesizedData.eccDegs;
   for k = 1:numel(eccRange)
       [~, idx(k)] = min(abs(eccDegs-eccRange(k)));
   end

   xRange = 1;
   eccDegs = eccDegs(idx);
   centerRadii = obj.synthesizedData.centerRadii(idx);
   surroundRadii = obj.synthesizedData.surroundRadii(idx); 
   centerPeakSensitivities = obj.synthesizedData.centerPeakSensitivities(idx);
   surroundPeakSensitivities = obj.synthesizedData.surroundPeakSensitivities(idx);
    
   hFig4 = figure(4); clf;
   theAxesGrid = plotlab.axesGrid(hFig4, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.02, ...
            'heightMargin', 0.02, ...
            'bottomMargin', 0.03, ...
            'topMargin', 0.01);
        
   idx = 0;
   for iRow = 1:2
       for iCol = 1:2
           theAxes = theAxesGrid{iRow,iCol};
           idx = idx+1;
           if (idx == 3)
               labelAxes = true;
           else
               labelAxes = false;
           end
           renderRF2DOutline(theAxes, eccDegs(idx), centerRadii(idx), surroundRadii(idx), xRange, labelAxes);
           drawnow;
       end
   end
   
   
   hFig5 = figure(5); clf;
   theAxesGrid = plotlab.axesGrid(hFig5, ...
            'rowsNum', 2, ...
            'colsNum', 2, ...
            'leftMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.02, ...
            'heightMargin', 0.02, ...
            'bottomMargin', 0.03, ...
            'topMargin', 0.01);
        
   peakSensitivityRange = [-200 1000];
   idx = 0;
   for iRow = 1:2
       for iCol = 1:2
           theAxes = theAxesGrid{iRow,iCol};
           idx = idx+1;
           if (idx == 3)
               labelAxes = true;
           else
               labelAxes = false;
           end
           renderRF1DOutline(theAxes, eccDegs(idx), centerRadii(idx), surroundRadii(idx), ...
               centerPeakSensitivities(idx), surroundPeakSensitivities(idx), xRange, peakSensitivityRange, labelAxes);
       end
   end
   
end

function renderRF1DOutline(theAxes,eccDegs, centerRadius, surroundRadius, centerPeakSensitivity, surroundPeakSensitivity, xRange, yRange, labelAxes)
    xs = surroundRadius*2;
    xs = linspace(-xs,xs,501);
    xc = centerRadius*2;
    xc = linspace(-xc,xc,501);
    c = centerPeakSensitivity * exp(-(xc/centerRadius).^2);  
    s =-surroundPeakSensitivity* exp(-(xs/surroundRadius).^2);

    area(theAxes,xc,c,'EdgeColor', 'r'); hold(theAxes, 'on');
    area(theAxes,xs,s, 'EdgeColor', 'b');
    axis(theAxes, 'square');
    set(theAxes, 'XLim', xRange*[-1 1], 'YLim', yRange, ...
       'XTick', [-1:0.5:1], 'YTick', [0], 'XColor', 'none', 'YColor', 'none');
    if (~labelAxes)
        set(theAxes, 'XTickLabel', {}, 'YTickLabel', {});
    end
    grid(theAxes, 'off');
    box(theAxes, 'off');
    text(theAxes, -1, yRange(2)*0.8, sprintf('%2.2f degs', eccDegs), 'FontSize', 17);
end

function renderRF2DOutline(theAxes,eccDegs, centerRadius, surroundRadius, xRange, labelAxes)
    angles = 0:5:360;
    xc = cosd(angles)*surroundRadius*2;
    yc = sind(angles)*surroundRadius*2;
    patch(theAxes, xc,  yc, [0.5 0.5 1]); hold(theAxes, 'on');
    xc = cosd(angles)*centerRadius*2;
    yc = sind(angles)*centerRadius*2;
    patch(theAxes, xc,  yc, [1 0.5 0.5]);
    axis(theAxes, 'square');
    set(theAxes, 'XLim', xRange*[-1 1], 'YLim', xRange*[-1 1], ...
       'XTick', [-1:0.5:1], 'YTick', [-1:0.5:1]);
    if (~labelAxes)
        set(theAxes, 'XTickLabel', {}, 'YTickLabel', {});
    end
    title(theAxes, sprintf('%2.2f degs', eccDegs));
end
