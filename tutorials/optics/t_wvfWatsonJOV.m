% Replicate the published curves in
%
%   Computing human optical point spread functions
%   A.B. Watson (2015), JOV
%
% See also
%


%%
ieInit;
callStack = dbstack;
isRunningTutorialTest = any(strcmp({callStack.name}, 'ieRunTutorialExampleTests'));

%%
d = displayCreate;
font = fontCreate('Z','Georgia',14,96);
scene = sceneCreate('letter', font, d);
scene = sceneSet(scene,'fov',0.5);
if ~isRunningTutorialTest
    sceneWindow(scene);
end

%%
% These are the indices and Values that Beau put in his first example.
%  (There are not equation numbers or page numbers, sigh).
%
%    psf = ZernikePointSpread[zc]
%
zIndices = [2,-2; 2,0; 2,2; ...
    3,-3; 3,-1; 3, 1; 3, 3; 
    4, -4; 4,-2; 4, 0; 4, 2; 4, 4];

zValues = [-0.0946, 0.0969, 0.305, ...
    0.0459, -0.121, 0.0264, -0.113, ...
    0.0292, 0.03, 0.0294, 0.0163, 0.064];


% Empty wavefront struct
wvf = wvfCreate;

for ii = 1:numel(zValues)
    jIndex = wvfZernikeNMToOSAIndex(zIndices(ii,1),zIndices(ii,2));
    wvf = wvfSet(wvf,'zcoeffs',zValues(ii),jIndex);    
end

% Save the original zcoeffs
zcoeffs = wvfGet(wvf,'zcoeffs');

%%  Not the same as his published figure
% Scale the errors up

if isRunningTutorialTest
    scaleFactors = 1;
else
    scaleFactors = [1, 2, 5];
end
for sFactor = scaleFactors
    wvf = wvfSet(wvf,'zcoeffs',sFactor*zcoeffs);
    wvf = wvfSet(wvf,'lcaMethod','human');
    wvf = wvfCompute(wvf);
    if ~isRunningTutorialTest
        wvfPlot(wvf,'psf','unit','um','wave',550,'plotrange',20);
    end
    if ~isRunningTutorialTest
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
    end
    oi = wvf2oi(wvf,'humanlens',true);
    oi.optics.lens.density = 0;
    oi = oiCompute(oi,scene,'pad value','mean');
    oi = oiSet(oi,'name',sprintf('Watson x %d',sFactor));
    if ~isRunningTutorialTest
        oiWindow(oi);
        oiPlot(oi,'psf',550); set(gca,'xlim',xlim,'ylim',ylim);
    end
end

% There is a numerical mis-match for DHB and me to figure out.
% The Z coefficients in Beau's paper and our calculations are off by
% about a factor of 5 somehow.  This is based on the amount of blur in
% the ImagePlot[blurredletter] figure with respect to arc min ('E')
if ~isRunningTutorialTest
    wvfPlot(wvf,'psf','unit','min','wave',550,'plotrange',10);
    wvfPlot(wvf,'image psf','unit','min','wave',550,'plotrange',10);
end

%% Diffraction limited approximation to eye
oi = oiCreate('diffraction limited');
oi = oiSet(oi,'optics fnumber', 5.67);
oi = oiSet(oi,'optics focal length',0.016);
oi = oiCompute(oi,scene,'pad value','mean');
name = oiGet(oi,'name'); oi = oiSet(oi,'name',sprintf('Diff %s',name));
if ~isRunningTutorialTest
    oiWindow(oi);
end

%% Marimont Wandell eye model
oi = oiCreate('human mw');
oi = oiCompute(oi,scene,'pad value','mean');
name = oiGet(oi,'name'); oi = oiSet(oi,'name',sprintf('MW %s',name));
if ~isRunningTutorialTest
    oiWindow(oi);
end

%% Thibos standard human
%
% Sharper than the Marimont/Wandell version
oi = oiCreate('wvf human');

% Needs a better interface. The lens transmittance is set to 1, but
% still we are getting blur due to chromatic aberration.
oi.optics.lens.density = 0;

oi = oiCompute(oi,scene,'pad value','mean');
name = oiGet(oi,'name'); oi = oiSet(oi,'name',sprintf('Thibos %s',name));
if ~isRunningTutorialTest
    oiWindow(oi);
end

if ~isRunningTutorialTest
    wvfPlot(wvf,'image wavefront aberrations','unit','um','wave',550);
end
zcoeffs = wvfGet(wvf,'zcoeffs');

%% Notice that this does not match aberration value in Watson

% The overall pattern and indices of the aberrations match.
%
% But it seems like we have a larger aberration using Beau's numbers.
% He shows a value of 1 or 2 microns, which seems right, and we have
% these huge values.
%
% If we scale the numerical value he uses to 0.63, instead of 63 given
% in the paper, we are much closer.  So, about a factor of 100
% difference somewhere in someone's units.
wvf = wvfCreate;
jIndex = wvfZernikeNMToOSAIndex(2,-2);
wvf = wvfSet(wvf,'zcoeff',0.63,jIndex);
if ~isRunningTutorialTest
    disp(wvfGet(wvf,'zcoeff'));
end
wvf = wvfSet(wvf,'lcaMethod','human');
wvf = wvfCompute(wvf);
if ~isRunningTutorialTest
    wvfPlot(wvf,'image wavefront aberrations','unit','um','wave',550)
    colormap("gray");
end

%% Show the wavefront aberrations for each Z coefficient up to 8s
if isRunningTutorialTest
    zernikeIndices = 2:3;
else
    zernikeIndices = 2:8;
end
for ii=zernikeIndices
    wvf = wvfCreate;
    wvf = wvfSet(wvf,'zcoeff',1,ii);
    wvf = wvfSet(wvf,'lcaMethod','human');
    wvf = wvfCompute(wvf);
    [n,m] = wvfOSAIndexToZernikeNM(ii);
    if ~isRunningTutorialTest
        wvfPlot(wvf,'image wavefront aberrations','unit','um','wave',550);
        colormap("gray"); title(sprintf('Z_%d^%d',n,m));
    end
end
