%
% RGCMosaicAnalyzer.visualize.coneExcitationsVsPhotocurrentResponse
%

function coneExcitationsVsPhotocurrentResponse(...
	theConeType, ...
	temporalSupportSeconds, ...
	theSingleConeExcitationsResponse, ...
	temporalSupportPhotocurrent, ...
	theConePhotoCurrentDifferentialResponse, ...
    theConeBackgroundPhotoCurrent, ...
    theConeExcitationsSingleConePeriodic, ...
    temporalSupportSecondsPeriodic)

	switch (theConeType)
        case cMosaic.LCONE_ID
          theColor = [231 187 192]/255;
        case cMosaic.MCONE_ID
          theColor = [128 207 218]/255;
        case cMosaic.SCONE_ID
          theColor = [0 0 1];
    end

      m1 = max(theSingleConeExcitationsResponse(:));
      m2 = min(theSingleConeExcitationsResponse(:));
      coneContrast = (m1-m2)/(m1+m2);
      dt = (temporalSupportSeconds(2)-temporalSupportSeconds(1))/2;

      coneExcitationsRange = [0 1800];
      coneExcitationsTicks = 0:100:2000;

      coneExcitationsRange = [0 500];
      coneExcitationsTicks = 0:50:2000;

      coneExcitationRateRange = [0 20000];
      coneExcitationRateTicks = 0:2000:30000;

      coneExcitationRateTickLabels = cell(1, numel(coneExcitationRateTicks));
      for i = 1:numel(coneExcitationRateTicks)
        coneExcitationRateTickLabels{i} = sprintf('%2.0fk', coneExcitationRateTicks(i)/1000);
      end

      photocurrentRange = [-75 -15];
      deltaPhotocurrentRange = [-25 25];

      %photocurrentRange = [-80 -0];
      %deltaPhotocurrentRange = [-45 35]

      photocurrentTicks = -80:5:80;
      timeTicks = 0:100:10000;

      dTSeconds = temporalSupportSeconds(2)-temporalSupportSeconds(1);

      hFig = figure(33); clf;
      set(hFig, 'Position', [10 10 1600 1150], 'Color', [1 1 1]);

      % The cone excitations response (1 period)
      subplot(2,3,3);

      stairs(temporalSupportSeconds*1e3, theSingleConeExcitationsResponse/dTSeconds, 'Color', 'k', 'LineWidth', 3.0);
      hold on;
      stairs(temporalSupportSeconds*1e3, theSingleConeExcitationsResponse/dTSeconds, 'Color', theColor, 'LineWidth', 1.5);
      plot([temporalSupportSeconds(1) temporalSupportSeconds(end)]*1e3, 0.5*(m1+m2)/dTSeconds*[1 1], '--', 'LineWidth', 1.0, 'Color', 'k', 'LineWidth', 1.0);
      ylabel('cone excitation rate (R*/sec)')
      xlabel('time (msec)');
      xtickangle(90);

      title(sprintf('response modulation: %2.2f', coneContrast), 'FontWeight', 'Normal');
      set(gca, 'YLim', coneExcitationRateRange, 'YTick', coneExcitationRateTicks, 'YTickLabel', coneExcitationRateTickLabels, ...
        'XTick', timeTicks, 'XLim', [temporalSupportSeconds(1) temporalSupportSeconds(end)+dt]*1e3, 'FontSize', 20);
      grid on; box off

      % The periodic cone excitations response
      subplot(2,3,[1 2]);

      stairs(temporalSupportSecondsPeriodic*1e3, theConeExcitationsSingleConePeriodic/dTSeconds, 'Color', 'k', 'LineWidth', 3);
      hold on;
      stairs(temporalSupportSecondsPeriodic*1e3, theConeExcitationsSingleConePeriodic/dTSeconds, 'Color', theColor, 'LineWidth', 1.5);
      plot([temporalSupportSecondsPeriodic(1) temporalSupportSecondsPeriodic(end)]*1e3, 0.5*(m1+m2)/dTSeconds*[1 1], '--', 'LineWidth', 1.0, 'Color', 'k', 'LineWidth', 1.0);
      ylabel('cone excitation rate (R*/sec)')
      xlabel('time (msec)');
      title(sprintf('response modulation: %2.2f', coneContrast), 'FontWeight', 'Normal');
      xtickangle(90);
      set(gca, 'YLim', coneExcitationRateRange, 'YTick', coneExcitationRateTicks, 'YTickLabel', coneExcitationRateTickLabels, ...
        'XTick', timeTicks, 'XLim', [temporalSupportSecondsPeriodic(1) temporalSupportSecondsPeriodic(end)]*1e3, 'FontSize', 20);
      grid on; box off


      % The photocurrent response (1 period)
      subplot(2,3,6);
      plot(temporalSupportPhotocurrent*1e3, theConePhotoCurrentDifferentialResponse, '-', 'Color', 'k', 'LineWidth', 3);
      hold on
      plot(temporalSupportPhotocurrent*1e3, theConePhotoCurrentDifferentialResponse, '-', 'Color', theColor, 'LineWidth', 1.5);
      plot(temporalSupportPhotocurrent*1e3,  temporalSupportPhotocurrent*0, 'k--');
      title(sprintf('(+): %2.2f, (-): %2.2f', max(theConePhotoCurrentDifferentialResponse), min(theConePhotoCurrentDifferentialResponse)), 'FontWeight', 'Normal');
      ylabel('differential pCurrent (pAmps)')
      xlabel('time (msec)');
      xtickangle(90);
      grid on; box off
      set(gca, 'YLim', deltaPhotocurrentRange, 'YTick', photocurrentTicks, 'XTick', timeTicks, 'XLim', [temporalSupportPhotocurrent(1) temporalSupportPhotocurrent(end)+dt]*1e3, 'FontSize', 20);

      subplot(2,3,[4 5]);
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentResponsePeriodic+theConeBackgroundPhotoCurrent, '-', 'Color', 'k', 'LineWidth', 3);
      hold on;
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentResponsePeriodic+theConeBackgroundPhotoCurrent, '-', 'Color', theColor, 'LineWidth', 1.5);
      plot(photocurrentResponseTimeAxisPeriodic*1e3, thePcurrentBackgroundResponseTransient, '--', 'Color', 'k', 'LineWidth', 1.0);
      set(gca, 'YLim', photocurrentRange , 'YTick', photocurrentTicks,  'XTick', timeTicks, 'XLim', [photocurrentResponseTimeAxisPeriodic(1) photocurrentResponseTimeAxisPeriodic(end)]*1e3, 'FontSize', 20);
      grid on; box off
      xtickangle(90);
      title(sprintf('(+): %2.2f, (-): %2.2f', max(theConePhotoCurrentDifferentialResponse), min(theConePhotoCurrentDifferentialResponse)), 'FontWeight', 'Normal');
      ylabel('pCurrent (pAmps)')
      xlabel('time (msec)');

end
