function varargout = v_oi(varargin)
%
% Test optical image creating functions
%
% Implicitly tests the opticsCreate functions, as well.
%
% Copyright Imageval LLC, 2009

% History:
%    08/31/23  dhb  This was passing but storing full structures.  I
%                   changed to do more computes and save the photons.  This will
%                   generalize better

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
    
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Initialize ISETBIO
    ieInit;

    % Create a scene to check oi function
    scene = sceneCreate;
    UnitTest.validationData('theScenePhotons',sceneGet(scene,'photons'));

    %% Diffraction limited simulation properties
    oi = oiCreate('diffraction limited');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'otf',[],550);
        oiPlot(oi,'otf',[],450);
    end
    % UnitTest.validationData('diffractionOI', oi);
    UnitTest.validationData('diffractionLimitedFromScenePhotons', oiGet(oi,'photons'));
    
    %% Wavefront (Thibos) human optics
    oi = oiCreate('wvf human');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    % UnitTest.validationData('humanWVF', oi);
    UnitTest.validationData('humanWVFFromScenePhotons', oiGet(oi,'photons'));

    %% A simple case used for testing
    oi = oiCreate('uniform ee');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    % UnitTest.validationData('EEoi', oi);
    UnitTest.validationData('unifromEEFromScenePhotons', oiGet(oi,'photons'));

    %% Make a scene and show some oiGets and oiCompute work
    oi = oiCreate('human');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'illuminance mesh linear');
    end
    %UnitTest.validationData('theScene',scene);
    %UnitTest.validationData('humanOIFromScene', oi);
    UnitTest.validationData('humanMWOIFromScenePhotons', oiGet(oi,'photons'));

    %% Check GUI control
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(oi);
        oiWindow;
        oiSet([],'gamma',1);
    end

end