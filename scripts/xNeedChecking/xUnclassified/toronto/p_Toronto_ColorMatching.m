%% p_Toronto_ColorMatching
%
% Illustrate some aspects of color matching through simulation.
%
% The script was used for a talk at the OSA.
%  * The first section creates a sample cone mosaic.
%  * The second section creates optics that match the human cornea/lens.
%  * The third section sets a peak radiance (you can choose the photon levels)
%    and wavelengths, and then produces cone absorptions in the human
%    sensors. 
%
% (c) ImagEval 2011

%% 
s_initISET
talkD = pwd;
saveFlag = 0;  % Don't save results

%%  Create a sample cone mosaic
params.sz = [128,192];
params.rgbDensities = [0.1 .6 .2 .1]; % Empty, L,M,S
params.coneAperture = [3 3]*1e-6;     % In meters
pixel = [];

% Here is a look at a typical cone mosaic.  Not used later.
cSensor = sensorCreate('human',pixel,params);
sensorConePlot(cSensor)

%%  Cone absorption as a function of mean signal level 

% Illustrate the separation in the LMS clouds for pairs of wavelengths
oi = oiCreate('human');

%%
pSize = 2e-6;   %Size in meters of the cone sampling apertures

% Create uniform fields on a human sensor (ideal) with Bayer sampling
% arrangement.  Then get the photon absorptions and plot them as a 3D graph
% in the next cell.
wSamples = [520  530];
% wSamples = 450:10:640;

% Number of row/col samples in the scene 
sz = 128;  

% It is possible to run this for a series of levels.  These levels are
% chosen to match typical levels of photopic illumination.
% It is easiest to run just one level, as per below.
%
% peakRadiance = [1e16 2e16 5e16 1e17 2e17 5e17 1e18];
peakRadiance = [ 5e16 5e17 ]*(30/100);

% We will make a series of scenes at different wavelengths and peak
% readiances. We will compute the sensor response.
scene  = cell(1,length(wSamples));
sensor = cell(1,length(wSamples));
for rr = 1:length(peakRadiance)
    for ww=1:length(wSamples)

        %% Create a monochromatic scene and set the radiance
        % The wavelength is specified in wSamples.
        scene{ww} = sceneCreate('uniform monochromatic',wSamples(ww),sz);
        scene{ww} = sceneSet(scene{ww},'peak photon radiance',peakRadiance(rr));
        % Some people scale for luminance, or equal L+M, which we could do.
        % We could go equal energy, not equal photon.
        % scene{ww} = sceneAdjustLuminance(scene{ww},100);

        % Compute the irradiance at the retina
        oi = oiCompute(scene{ww},oi);
        sceneGet(scene{ww},'mean luminance')

        % Create a human sensor.
        sensor{ww} = sensorCreate('ideal',[],pSize,'human','bayer');
        % Integrate for 100 ms.
        sensor{ww} = sensorSet(sensor{ww},'exposure time',0.10);
        
        % Compute the sensor absorptions
        sensor{ww} = sensorCompute(sensor{ww},oi);
        
        % Give the scene a name
        sensor{ww} = sensorSet(sensor{ww},'name',sprintf('wave %.0f',wSamples(ww)));
        
        % If you want to have a look at the image, run this line.
        % vcAddAndSelectObject(sensor{ww}); sensorImageWindow;
    end

    %% Extract the data for plotting
    L = cell(1,length(wSamples));
    M = cell(1,length(wSamples));
    S = cell(1,length(wSamples));

    % The cones are in slots 2-4 when not ideal
    slot = [ 1 2 3];  % Ideal case
    % slot = [2 3 4];   % Typical human 1621 case.  1 empty, 6 L, 2 M, 1 S
    for ww=1:length(wSamples)
        L{ww} = sensorGet(sensor{ww},'electrons',slot(1));
        M{ww} = sensorGet(sensor{ww},'electrons',slot(2));
        S{ww} = sensorGet(sensor{ww},'electrons',slot(3));
        n = min(100,length(L{ww}));

        % Make the same length
        S{ww} = S{ww}(1:n); M{ww} = M{ww}(1:n); L{ww} = L{ww}(1:n);
    end

    % Plot the absorptions
    f = vcNewGraphWin;
    % sym = {'b.','g.','r.','c.','bs','gs','rs','cs','bo','go','ro','co','bx','gx','rx','cx'};
    sym = {'b.','g.','r.','c.','k.'};
    az = 65.5; el = 30;
    for ww=1:length(wSamples)
        s = mod(ww,length(sym))+1;
        plot3(L{ww}(:),M{ww}(:),S{ww}(:),sym{s})
        view([az el])
        hold on
    end

    xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on

    % Probably don't want to save in general
    if saveFlag
        level = 10*round(log10(peakRadiance(rr))/10);
        fname = fullfile(talkD,'images','CMF',sprintf('CMF-%.0f-%.1f.fig',wSamples(1),level));
        saveas(f,fname,'fig');
    end
end

% vcAddAndSelectObject(sensor{ww}); sensorImageWindow;

