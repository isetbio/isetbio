function varargout = v_oi(varargin)
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
    s_initISET;

    %% Diffraction limited simulation properties
    oi = oiCreate;
    if (runTimeParams.generatePlots)
        oiPlot(oi,'otf',[],550);
        oiPlot(oi,'otf',[],450);
    end
    UnitTest.validationData('diffractionOI', oi);

    %% Human optics
    oi = oiCreate('human');
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    UnitTest.validationData('humanOI', oi);

    %% Make a scene and show some oiGets and oiCompute work
    scene = sceneCreate;
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'illuminance mesh linear');
    end
    UnitTest.validationData('theScene',scene);
    UnitTest.validationData('humanOIFromScene', oi);


    %% Check GUI control
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(oi);
        oiWindow;
        oiSet([],'gamma',1);
    end

    %% End
end