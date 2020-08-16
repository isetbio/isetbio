function plotDeconvolutionModel(deconvolutionModel)
    
    visualizeCenterDeconvolutionModel(deconvolutionModel.center);
    visualizeSurroundDeconvolutionModel(deconvolutionModel.surround);
    
end

function visualizeSurroundDeconvolutionModel(deconvolutionModel)
    hFig = figure(2);
    set(hFig, 'Name', sprintf('RF surround deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    rowsNum = 2;
    colsNum = numel(deconvolutionModel.eccDegs);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum);
        
    maxSurroundPeakSensitivity = 1; % max([max(deconvolutionModel.visualGain(:)) max(deconvolutionModel.retinalGain(:))]);
    maxSurroundCharacteristicRadius = max([max(deconvolutionModel.visualCharacteristicRadius(:)) max(deconvolutionModel.retinalCharacteristicRadius(:))]);
    
    for eccIndex = 1:numel(deconvolutionModel.eccDegs)
        % Visualize how we estimate the retinal characteristic radius of
        % the surround from the visual characteristic radius of the surround
        ax = theAxesGrid{1, eccIndex};
        scatter(ax, deconvolutionModel.visualCharacteristicRadius(eccIndex,:), ...
                    deconvolutionModel.retinalCharacteristicRadius(eccIndex,:), ...
                    'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        hold(ax, 'on');
        plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1.0);
        set(ax, 'XLim', [0 maxSurroundCharacteristicRadius], 'YLim', [0 maxSurroundCharacteristicRadius]);
        axis(ax, 'square');
        if (eccIndex == 1)
            ylabel(ax, 'surround retinal characteristic radius');
        else
            set(ax, 'YTickLabel', {})
        end
        xlabel(ax, 'surround visual characteristic radius');
        title(ax,sprintf('ecc: %2.1f degs', deconvolutionModel.eccDegs(eccIndex)));
          
        % Visualize how we estimate the retinal peak sensitivity of the surround 
        % from the visual peak sensitivity of the surround
        ax = theAxesGrid{2, eccIndex};
        scatter(ax, deconvolutionModel.visualGain(eccIndex,:), ...
                    deconvolutionModel.retinalGain(eccIndex,:), ...
                    'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        hold(ax, 'on');
        plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1.0);
        
        set(ax, 'XLim', [0 maxSurroundPeakSensitivity], 'YLim', [0 maxSurroundPeakSensitivity], ...
            'XTick', 0:0.1:1, 'YTick', 0:0.1:1);
        axis(ax, 'square');
        if (eccIndex == 1)
            ylabel(ax, 'surround retinal peak sensitivity');
        else
            set(ax, 'YTickLabel', {})
        end
        xlabel(ax, 'surround visual peak sensitivity');
    end

        
end

function visualizeCenterDeconvolutionModel(deconvolutionModel)

    w = WatsonRGCModel();
    retinalNasalEccentricitiesDegs = logspace(log10(0.1), log10(30), 100);
    coneRFSpacingsDegs = w.coneRFSpacingAndDensityAlongMeridian(retinalNasalEccentricitiesDegs, ...
            'nasal meridian','deg', 'deg^2', ...
            'correctForMismatchInFovealConeDensityBetweenWatsonAndISETBio', false);
    coneApertureDiameterDegs = WatsonRGCModel.coneApertureToDiameterRatio * coneRFSpacingsDegs;
    coneRadiiDegs = 0.5*coneApertureDiameterDegs;
    
    ck = CronerKaplanRGCModel('generateAllFigures', false);
    ck.setupPlotLab(0, 25, 15);
    
    hFig = figure(1); clf;
    set(hFig, 'Name', sprintf('RF center deconvolution model for subject %d and ''%s'' quadrant.', ...
        deconvolutionModel.subjectID, deconvolutionModel.quadrant));
    
    rowsNum = 2;
    colsNum = numel(deconvolutionModel.eccDegs);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.01, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum);
        
    maxCenterPeakSensitivity = 1; %max([max(deconvolutionModel.visualGain(:)) max(deconvolutionModel.retinalGain(:))]);
  
    for eccIndex = 1:numel(deconvolutionModel.eccDegs)
        [~,idx]  = min(abs(retinalNasalEccentricitiesDegs-deconvolutionModel.eccDegs(eccIndex)));
        coneRadiusDegsAtClosestEccentricity = coneRadiiDegs(idx);
        
        % Visualualize how we estimate the visual characteristic radius of the RF center from the # of cone in RF center
        ax = theAxesGrid{1, eccIndex};
        
        scatter(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualCharacteristicRadiusMin(eccIndex,:), ...
            'r^', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        hold(ax,'on');
        
        scatter(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualCharacteristicRadiusMax(eccIndex,:), ...
            'rv', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        
        % The mean (min+max)/2 characteristic radius of the VISUAL center
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), ...
              0.5*(deconvolutionModel.visualCharacteristicRadiusMax(eccIndex,:)+deconvolutionModel.visualCharacteristicRadiusMin(eccIndex,:)), ...
              'r-', 'LineWidth', 2);
        
        % The mean (min+max)/2 characteristic radius of the RETINAL center
        plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), ...
            0.5*(deconvolutionModel.retinalCharacteristicRadiusMin(eccIndex,:) +  deconvolutionModel.retinalCharacteristicRadiusMax(eccIndex,:)), ...
              'b--', 'LineWidth', 1.5);
        
        % The characteristic radius of a single cone
        plot(ax,retinalNasalEccentricitiesDegs, retinalNasalEccentricitiesDegs*0 + coneRadiusDegsAtClosestEccentricity, ...
            'k--', 'LineWidth', 1.5);
         
        set(ax, 'XLim', [0.9 31], 'YLim', [0.001 1], ...
            'YTick', [0.003 0.01 0.03 0.1 0.3 1], 'YTickLabel', {'0.003', '0.01', '0.03', '0.1', '0.3', '1.0'}, ...
            'XTick', [1 2 3 4 6 10 15 30 100], ...
            'XScale', 'log', 'YScale', 'log');
        axis(ax, 'square');
        
        if (eccIndex == 1)
            ylabel(ax, ['center characteristic radius ({\color{red} visual, \color{blue} retinal})']);
        else
            set(ax, 'YTickLabel', {})
        end
        xlabel(ax, '# of cones in RF center');
         
        title(ax,sprintf('ecc: %2.1f degs', deconvolutionModel.eccDegs(eccIndex))); 
        
        
        % Visualize how we estimate retinal peak sensitivity of the RF
        % center from the visual peak sensitivity of the RF center
        ax = theAxesGrid{2, eccIndex};
        scatter(ax, deconvolutionModel.visualGain(eccIndex,:), ...
                    deconvolutionModel.retinalGain(eccIndex,:), ...
                    'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [0.8 0.3 0.3]);
        hold(ax, 'on');
        plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1.0);
        set(ax, 'XLim', [0 maxCenterPeakSensitivity], 'YLim', [0 maxCenterPeakSensitivity], ...
             'XTick', 0:0.1:1, 'YTick', 0:0.1:1);
        axis(ax, 'square');
        if (eccIndex == 1)
            ylabel(ax, 'center retinal peak sensitivity');
        else
            set(ax, 'YTickLabel', {})
        end
        xlabel(ax, 'center visual peak sensitivity');
        
        
%         % Visual gain attenuation as a function of number of cone in RF center
%         ax = theAxesGrid{2, eccIndex};
%         plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.visualGain(eccIndex,:), ...
%             'r-', 'MarkerFaceColor', [1 0.5 0.5]);
%         hold(ax, 'on');
%         plot(ax,deconvolutionModel.centerConeInputsNum(eccIndex,:), deconvolutionModel.retinalGain(eccIndex,:), ...
%             'b-', 'LineWidth', 1.5);
%         set(ax, 'XLim', [0.9 31], 'YLim', [0 1.04], 'XScale', 'log', ...
%             'XTick', [1 2 3 4 6 10 15 30 100])
%         xlabel(ax,'# of center cones');
%         if (eccIndex == 1)
%             ylabel(ax, ['center peak gain ({\color{red} visual, \color{blue} retinal})']);
%         else
%             set(ax, 'YTickLabel', {})
%         end
        drawnow
   end
               
end