function bipolarFilt = bipolarFilter(obj, cmosaic, varargin)
% Estimate bipolar temporal filter 
%
%    bipolarFilt = bipolarFilter(obj, cmosaic, 'graph',logical)
% 
% This routine calculates the bipolar impulse response function such that
% the convolution of the outersegment and the bipolar equals the impulse
% response of the rgc.
%
% BW: Mysterious, what we are doing here.  With this calculation, we should
% not further apply the RGC impulse response.  I think.  To discuss with
% the group what the IR of the RGC means in the ir.mosaic{1}.tCenter field.
% 
% Example:
%   cMosaic = coneMosaic;
%   cMosaic.os.linearFilters(cMosaic);
%   bp = bipolarCreate(cMosaic);
%   bpFilter = bipolarFilter(bp, cMosaic,'graph',true);
%
%  If you used a cMosaic with a real time varying signal, then you can do
%  this.  See t_bpTemporal
%
%   vcNewGraphWin; plot(cMosaic.timeAxis,bpFilter,'o-');
%
% (c) isetbio team JRG/BW 12/2016

%% parse input parameters
p = inputParser;

p.addRequired('obj', @(x) (isa(x, 'bipolar')));
p.addRequired('cmosaic', @(x) (isa(x, 'coneMosaic')));  
p.addParameter('graph',false,@islogical);

p.parse(obj, cmosaic, varargin{:});
graph = p.Results.graph;

%% Get outersegment impulse response filters

% Get the finest time sampling from the os for the calculation.  We
% downsample in time at the end, using cmosaic.timeAxis
if isa(cmosaic.os,'osLinear')
    osFilt   = cmosaic.os.lmsConeFilter(:,1);
    osTime   = cmosaic.os.timeAxis;
else
    osFiltTemp = linearFilters(cmosaic.os, cmosaic);
    osFilt = osFiltTemp(:,1);
    osTime   = ((1:size(osFilt,1)) - 1) * cmosaic.os.timeStep;
end


%% Create the Pillow impulse response function that we match

% We match the impulse response of an inner retina type irPhys, the values
% for eyeSide, eyeRadius and eyeAngle have no effect, because those are not
% dependent on the properties of the retinal piece used in the Chichilnisky
% Lab experiment.
% clear params
% params.name = 'dummy';
% params.eyeSide = 'left'; 
% params.eyeRadius = 0; 
% params.eyeAngle  = 0; 
% cellType = 'onParasol';   % Could be offParasol, or ...
% 
% % First the inner retina and then a mosaic
% % The inner retina seems to use the obj (bipolar) time base to set up the
% % impulse response in the mosaicCreate.  Look into this and get it right,
% % sigh.
% innerRetina = ir(obj, params);
% innerRetina.mosaicCreate('type',cellType,'model','GLM');
% 
% % The temporal impulse response function of the mosaic.
% cellNum = 1;
% rgcFilt = innerRetina.mosaic{1}.tCenter{cellNum};
% rgcTime = innerRetina.mosaic{1}.timeAxis;

% I chose these parameters so that the bipolar IR would be delayed past the
% OS impulse response.
params.filterDuration = 0.4; 
params.samplingTime = 0.001;
params.cellType = obj.cellType;
[rgcFilt,rgcTime ] = rgcImpulseResponsePillow(params);
% vcNewGraphWin; plot(rgcTime,rgcFilt);

% Put the two IR functions on the same very fine time base
timeBase = (0:4095)*1e-4;    % Time samples in seconds.  About 400 ms total
rgcFilt = interp1(rgcTime,rgcFilt,timeBase,'linear',0);
osFilt  = interp1(osTime, osFilt,timeBase,'linear',0);


%% Find the bipolar impulse response

% Convolve the estimated bipolar IR with a Gaussian window for low-pass
% filter to reduce the noise.  The size of the gaussian window controls the
% temporal frequencies that we retain.  The circshift of the Gaussian
% controls the onset of the impulse response.  We shift by an amount the
% keeps the filter causal.
gaussVar = 100;
gw1 = gausswin(length(rgcFilt),gaussVar)';
gw1 = circshift(gw1,-round(length(rgcFilt)/2) + 0*gaussVar,2);
gw1 = gw1/sum(gw1);   % Unit area for no DC amplification 

% Strange indeterministic error here - occasionally, depending on the
% average isomerizations computed in @osLinear/linearFilters.m, the
% numerator and denominator of fftBipolarFilt are 0 at the same point,
% which results in NaN. Here we add a small epsilon 1e-15 to denominator if
% this happens in order to get rid of the NaN values.
fftBipolarFilt = fft(gw1) .* fft((rgcFilt)) ./ (fft((osFilt)));
nanInd = isnan(fftBipolarFilt); 
if sum(nanInd(:))>0
    fftBipolarFiltEPS = fft(gw1) .* fft((rgcFilt)) ./ (1e-15+fft((osFilt)));
    fftBipolarFilt((nanInd)) = fftBipolarFiltEPS((nanInd));
end
bipolarFilt = ifft(fftBipolarFilt);

% Graph the various curves (except the Gaussian)
if graph
    vcNewGraphWin([],'wide'); 
    subplot(1,2,1)
    plot(timeBase,osFilt/max(osFilt(:)),'b-',...
        timeBase,rgcFilt/max(rgcFilt(:)),'r-',...
        timeBase,bipolarFilt/max(bipolarFilt(:)),'g-','linewidth',1);
    grid on; set(gca,'xlim',[-0.05 0.4]); xlabel('Time (sec)');
    set(gca,'fontsize',14); legend({'os','rgc','bipolar'})
    title('Normalized peak filters')
    
    % Estimate the rgcFilt from the osFilt and bipolarFilt
    subplot(1,2,2)
    est = ifft(fft(osFilt) .* fft(bipolarFilt));
    plot(timeBase,est,'r-',timeBase,rgcFilt,'k--');
    legend('estimate','rgcFilt'); grid on;
    set(gca,'xlim',[-0.05 0.4]); xlabel('Time (sec)');
end


switch obj.cellType
    case{'ondiffuse','onmidget','onsbc','sbc'}
        signVal = 1;
    otherwise
        signVal = -1;
end
            
% Down sample the time axis to match the cone mosaic sample times
bipolarFilt = signVal*interp1(timeBase,bipolarFilt,cmosaic.timeAxis,'linear',0);

end