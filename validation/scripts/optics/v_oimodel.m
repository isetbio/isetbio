function varargout = v_oimodel(varargin)
%
% Test optical image functions
%
%
% Copyright Imageval LLC, 2009

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Initialize ISETBIO
    ieInit;

    %% Diffraction limited simulation properties
    oi = oiCreate('diffraction limited');
    if (runTimeParams.generatePlots)
        oiPlot(oi,'otf',[],550);
        oiPlot(oi,'otf',[],450);
    end
    UnitTest.validationData('diffractionOI', oi);

    %% Shift invariant 
    oi = oiCreate('shift invariant');
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    UnitTest.validationData('humanOI', oi);

    %% Check GUI control
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(oi);
        oiWindow;
        oiSet([],'gamma',1);
    end

end