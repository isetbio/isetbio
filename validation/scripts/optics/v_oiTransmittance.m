function varargout = v_oiTransmittance(varargin)
%
% Validate some of the optical image transmittance calculations

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
%ieInit;

% Tolerance fraction
toleranceFraction = 0.005;

%% Human optics
oi = oiCreate;
transmittance1 = oiGet(oi,'optics transmittance');
wavelenghts1 = oiGet(oi,'optics transmittance wave');
theTolerance = mean(transmittance1(:))*toleranceFraction;
UnitTest.validationData('wavelengths1',wavelenghts1);
UnitTest.validationData('transmittance1',transmittance1, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'transmittance1', theTolerance);

if (runTimeParams.generatePlots)
    vcNewGraphWin;
    plot(wavelenghts1,transmittance1);
    xlabel('Wave'); ylabel('Transmittance'); grid on
end

%% Create a scene and check that its wavelength is imposed on transmittance
wavelengths2 = 450;
scene = sceneCreate;
scene = sceneSet(scene,'wavelength',wavelengths2);
oi = oiCompute(oi,scene);

% Notice that the new get respects the new computation
transmittance2 = oiGet(oi,'optics transmittance');
theTolerance = mean(transmittance2(:))*toleranceFraction;
UnitTest.validationData('wavelengths2',wavelengths2);
UnitTest.validationData('transmittance2',transmittance2, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'transmittance2', theTolerance);
if (runTimeParams.generatePlots)
    oiShowImage(oi);
end

%% Diffraction optics
%
% Notice that we can pull out transmittance at custom wavelengths
% that differ from those in the scene. 
scene = sceneCreate;
scene = sceneSet(scene,'wavelength',550);
oi = oiCreate('diffraction');
oi = oiCompute(oi,scene);
wavelengths3 = 400:10:500;
transmittance3 = oiGet(oi,'optics transmittance',wavelengths3);
theTolerance = mean(transmittance3(:))*toleranceFraction;
UnitTest.validationData('wavelengths3',wavelengths3);
UnitTest.validationData('transmittance3',transmittance3, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'transmittance3', theTolerance);
if (runTimeParams.generatePlots)
    vcNewGraphWin;
    plot(wavelenghts1,transmittance1);
    xlabel('Wave'); ylabel('Transmittance'); grid on
    oiShowImage(oi);
end

end
