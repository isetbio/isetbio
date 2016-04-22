% t_rgcSubunitHochShap
% 
% Demonstrates the inner retina object calculation for the subunit RGC
% model (from Hochstein & Shapley, 1976)
% 
% Nonlinear spatial subunits within RGC RF. After the spatial convolution,
% the signal is put through a rectifying linearity (at zero).
% 
% 3/2016 BW JRG HJ (c) isetbio team

%%
ieInit

%% Movie of the a monochromatic region of rod absorptions (TODO)


% For the moment, we take what's up there which is a bunch of cone
% isomerizations.  We will make a directory for creating rod stimuli as
% well.

% Get data from isetbio archiva server
rd = RdtClient('isetbio');
rd.crp('/resources/data/istim');
a = rd.listArtifacts;

% Pull out .mat data from artifact
whichA =1 ;
data = rd.readArtifact(a(whichA).artifactId);
% iStim stores the scene, oi and cone absorptions
iStim = data.iStim;
absorptions = iStim.absorptions;

% Grating subunit stimulus
% params.barWidth = 24;
% iStim = ieStimulusGratingSubunit;
% absorptions = iStim.absorptions;

% White noise
% iStim = ieStimulusWhiteNoise;

% Show raw stimulus for osIdentity
coneImageActivity(absorptions,'dFlag',true);

%% Photocurrent calculation

% We have the moment by moment absorptions.  We now want to create the
% tempmorally filtered version.  If G&M had given us a temporal impulse
% response for the photoreceptor, we would have used it.  For this
% calculation, we use the ISETBIO default.
os = osCreate('linear');

% Set size of retinal patch
patchSize = sensorGet(absorptions,'width','um');
os = osSet(os, 'patch size', patchSize);

% Set time step of simulation equal to absorptions
timeStep = sensorGet(absorptions,'time interval','sec');
os = osSet(os, 'time step', timeStep);

% Plot the photocurrent for a center pixel
%
%   osPlot(os,absorptions);
% 

% This are the Rieke et al. cone impulse response functions
%   osPlot(os,absorptions,'filters')
%   grid on
%


%%  EJ says for the bipolar model we can use this temporal impulse
%
%

% Get the filter from the osLinear

% Take the derivative

% add the derive (scaled a bit) back into the impulse response

% Put the whole thing back into the os object


% Set osI data to raw pixel intensities of stimulus
os = osCompute(os,absorptions);


%% G&M tell us to half-wave rectify the photocurrent response.

% What the hell does that mean?  Where is zero?
% We should write a hwRect function
%
%   out = hwRect(data,val);
%
eZero = -50;
hwrCurrent = ieHwrect(os.coneCurrentSignal,eZero);

%% Then we need a little spatial summation over the cones

% This is like a bipolar cell, but actually it could be the same code as in
% the spatial summation of the RGC
%
%  bipolar = spatialSummation(hwr,params);
%
%  spatialTemporalSummation()
%

kernel = fspecial('gaussian',[9,9],3);

bipolar = ieSpaceTimeFilter(hwrCurrent,kernel);

% For visualization, set the bipolar current to positive
% Not working correctly!  Try to understand how to visualize positive and
% negative numbers.  Maybe voltImageActivity ... that is thought through
% correctly for positive and negative numbers.
% bmosaic = sensorSet(absorptions,'photons',bipolar);
% coneImageActivity(bmosaic,'dFlag',true);

% Not sure if this is detailed enough; not used after this point.
strideSubsample = 4;
bipolarSubsample = ieImageSubsample(bipolar, strideSubsample);

%% Show bipolar activity
figure;
for frame1 = 1:size(bipolar,3)
    imagesc(squeeze(bipolar(:,:,frame1)));
    colormap gray; drawnow;
end
close;

%% Build the inner retina object
% Still working on bipolar processing
% Need to determine how many bipolars per RGC in order to subsample

clear params innerRetina0
params.name      = 'Macaque inner retina 1'; % This instance
params.eyeSide   = 'left';   % Which eye
params.eyeRadius = 4;        % Radius in mm
params.eyeAngle  = 90;       % Polar angle in degrees

innerRetina0 = irCreate(os, params);

% Create a coupled GLM model for the on midget ganglion cell parameters
innerRetina0.mosaicCreate('model','lnp','type','on midget');


%% Calculate bipolar inputs by subsampling full bipolar response

% Half-wave rectify
eZero = -50;
hwrBipolar = ieHwrect(bipolar,eZero);

bipolarSize = size(hwrBipolar);

% Assume 2x2 bipolars for each RF
% See Freeman, Field, Sher, Litke, Simoncelli, Chichilnisky, eLife 2016
% http://elifesciences.org/content/4/e05241v1

numberSubunits = [2 2];%mosaic.numberSubunits;

numberRGCs = size(innerRetina0.mosaic{1}.cellLocation);

bipolarRFsize = floor(bipolarSize(1:2)./numberRGCs);

% The size of each bipolar spatial input 
suSize1 = floor(bipolarRFsize(1)/numberSubunits(1));
suSize2 = floor(bipolarRFsize(2)/numberSubunits(2));

% Do subsampling
for rgcInd1 = 1:numberRGCs(1)
    for rgcInd2 = 1:numberRGCs(2)
        rgcBipolarSum{rgcInd1,rgcInd2} = 0;      
        
        for suInd1 = 1:numberSubunits(1)
            for suInd2 = 1:numberSubunits(2)               
                
                % Find center of subunit and sample value from there
                xCenter = (rgcInd1-1)*bipolarRFsize(1) + (suInd1)*round(suSize1/2);
                yCenter = (rgcInd2-1)*bipolarRFsize(2) + (suInd2)*round(suSize2/2);
                
                bipolarCenter{rgcInd1,rgcInd2}(suInd1,suInd2,:) = hwrBipolar(xCenter,yCenter,:);
                
                % Add to RGC input
                rgcBipolarSum{rgcInd1,rgcInd2} = rgcBipolarSum{rgcInd1,rgcInd2} + squeeze(bipolarCenter{rgcInd1,rgcInd2}(suInd1,suInd2,:))';
            end
        end
        
        
    end
end

% Set bipolar output to rgc linear input
mosaicSet(innerRetina0.mosaic{1},'responseLinear', rgcBipolarSum);
irPlot(innerRetina0,'linear');
%% Compute RGC mosaic responses

innerRetina0 = irCompute(innerRetina0, os);
irPlot(innerRetina0, 'psth');
irPlot(innerRetina0, 'linear');
% irPlot(innerRetina0, 'raster');

%% Show me the PSTH for one particular cell

% irPlot(innerRetina0, 'psth response','cell',[2 2]);
% irPlot(innerRetina0, 'raster','cell',[1 1]);