function generateMidgetRF(obj, cellID, visualizeDynamics)

    obj.spatialRF = struct();

    % Get temporal transfer params for the cellID
    params = obj.temporalKernelParams(cellID);

    % Generate temporal impulse response function for the center
    centerTemporalTransferFunction = obj.generateTemporalTrasferFunction(params.center);
    [obj.temporalRF.center.impulseResponse, obj.temporalRF.center.temporalSupport] = ...
        RGCRF.impulseResponseFromOneSidedTransferFunction(centerTemporalTransferFunction, obj.tfSupport);

    % Generate temporal impulse response function for the surround
    surroundTemporalTransferFunction = obj.generateTemporalTrasferFunction(params.surround);
    [obj.temporalRF.surround.impulseResponse, obj.temporalRF.surround.temporalSupport] = ...
        RGCRF.impulseResponseFromOneSidedTransferFunction(surroundTemporalTransferFunction, obj.tfSupport);


    if (visualizeDynamics)
        temporalTransferFunctionMagnitude = abs(centerTemporalTransferFunction);
        temporalTransferFunctionPhase = angle(centerTemporalTransferFunction);
        figure(1);
        ax = subplot(1,3,1);
        plot(ax,obj.tfSupport , temporalTransferFunctionMagnitude, 'r-', 'LineWidth', 1.5);
        set(ax, 'Yscale', 'log', 'XLim', [0.5 64], 'XScale', 'log', 'XTick', [0.25 0.5 1 2 4 8 16 32]);
        xlabel(ax,'temporal frequency (Hz)')
        grid(ax, 'on');
        
        ax = subplot(1,3,2);
        plot(ax,obj.tfSupport , unwrap(temporalTransferFunctionPhase)/pi, 'r-', 'LineWidth', 1.5);
        set(ax, 'XLim', [0.5 64], 'XScale', 'log', 'XTick', [0.25 0.5 1 2 4 8 16 32 64], 'YLim', [-6 1]);
        xlabel(ax,'temporal frequency (Hz)')
        ylabel(ax,'phase (pi radians)')
        grid(ax, 'on')

        temporalTransferFunctionMagnitude = abs(surroundTemporalTransferFunction);
        temporalTransferFunctionPhase = angle(surroundTemporalTransferFunction);
        ax = subplot(1,3,1);
        hold(ax, 'on');
        plot(ax, obj.tfSupport , temporalTransferFunctionMagnitude, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');

        ax = subplot(1,3,2);
        hold(ax, 'on');
        plot(ax, obj.tfSupport , unwrap(temporalTransferFunctionPhase)/pi, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');

        ax = subplot(1,3,3);
        plot(ax, obj.temporalRF.center.temporalSupport*1e3, obj.temporalRF.center.impulseResponse, 'r-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, obj.temporalRF.surround.temporalSupport*1e3, obj.temporalRF.surround.impulseResponse, 'b-', 'LineWidth', 1.5);
        hold(ax, 'off');
        set(ax, 'XLim', [0 200], 'XTick', 0:25:200);
        grid(ax, 'on');
        xlabel(ax,'time (msec)')
    end
end




