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
    % ieInit;

    % Tolerance fraction
    toleranceFraction = 0.005;

    % Create a scene to check oi function
    scene = sceneCreate;
    theScenePhotons = sceneGet(scene,'photons');
    theTolerance = mean(theScenePhotons(:))*toleranceFraction;
    UnitTest.validationData('theScenePhotons',theScenePhotons, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'theScenePhotons', theTolerance);

    %% Diffraction limited simulation properties
    oi = oiCreate('diffraction limited');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'otf',[],550);
        oiPlot(oi,'otf',[],450);
    end
    theOiPhotons = oiGet(oi,'photons');
    theTolerance = mean(theOiPhotons(:))*toleranceFraction;
    UnitTest.validationData('diffractionLimitedFromScenePhotons', theOiPhotons, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'diffractionLimitedFromScenePhotons', theTolerance);    

    %% Wavefront (Thibos) human optics
    oi = oiCreate('wvf human');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    theOiPhotons = oiGet(oi,'photons');
    theTolerance = mean(theOiPhotons(:))*toleranceFraction;
    UnitTest.validationData('humanWVFFromScenePhotons', theOiPhotons, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'humanWVFFromScenePhotons', theTolerance);    

    %% A simple case used for testing
    oi = oiCreate('uniform ee');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'psf',[],420);
        oiPlot(oi,'psf',[],550);
    end
    theOiPhotons = oiGet(oi,'photons');
    theTolerance = mean(theOiPhotons(:))*toleranceFraction;
    UnitTest.validationData('unifromEEFromScenePhotons', theOiPhotons, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'unifromEEFromScenePhotons', theTolerance);    

    %% Make a scene and show some oiGets and oiCompute work
    oi = oiCreate('human');
    oi = oiCompute(oi,scene);
    if (runTimeParams.generatePlots)
        oiPlot(oi,'illuminance mesh linear');
    end
    theOiPhotons = oiGet(oi,'photons');
    theTolerance = mean(theOiPhotons(:))*toleranceFraction;
    UnitTest.validationData('humanMWOIFromScenePhotons', theOiPhotons, ...
        'UsingTheFollowingVariableTolerancePairs', ...
        'humanMWOIFromScenePhotons', theTolerance);    

    %% Check GUI control
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(oi);
        oiWindow;
        oiSet([],'gamma',1);
    end

end