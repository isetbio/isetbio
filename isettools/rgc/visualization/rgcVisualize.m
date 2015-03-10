function rgcVisualize(plotType,rgcP,varargin)
% Visualization routines for internal variables
%
%   rgcVisualize(plotType,rgcP,varargin)
%
% Gateway routine for plotting rgcParameters/Layers
%
% The plotTypes are
%   ** RGC Properties **
%     'center tirf' - Temporal impulse response of the center
%     'surround tirf' - Off the surround
%     'all tirf'      - Sum, center and surround TIRFs, plotted together
%     'feed back tr' - Feedback temporal response
%     'coupling tr' - Coupling temporal response
%
%   ** Images **********
%    'RF Image'         - receptive field of single RGC
%    'RF Center'        - receptive field of single RGC, center only
%    'RF Surround'      - receptive field of single RGC, surround only
%    'linear TS mean'   - mean of linear time series over time
%    'spikes mean'      - mean spike rate (blurring can be set(??))
%    'spikes mean mesh' - mesh of the mean spike rate (blurring can be set (?))
%        To control the blurring, which is usually a support of 7 and spread of
%        3, you can use a vector for the 4th argument, say [ 5, 2].
%
%   ** Movies **********
%    'Stimulus movie'   - NYI. Read the scene from cones.scene.
%                         Read the eye positions from sensor.
%                         Make a movie of the scene with eye positions.
%    'RF movie'         - NYI. Series of spatially blurred cone images?
%    'Cone Voltages'    - movie of cone voltages
%
%    'Linear TS'        - Movie of Linear time series, including spatial blur, eye movement
%                           temporal impulse response
%    'RGC TS'           - Movie of RGC time series (prior to thresholding)
%    'Spikes TS'        - Movie of RGC spikes (binary)
%
%   ** Other **********
%    'connectm'  - Just a print out of the connections in the middle
%
% Consider: lints and lintsimage are somewhat redundant. we could replace
% all such pairs with a two argument call:
%    either, rgcVisualize('linTS',rgcP, 'movie')
%    or      rgcVisualize('linTS',rgcP, 'image')
% Etc.
%
% Examples:
%
%   MANY PROBLEMS HERE .... FIX UP
%
% ... if you have computed absorptions for an rgc class (rgcP) ...
%    [spikes, rgcData] = rgcComputeSpikes(absorptions,rgcP);
%    rgcVisualize('linTS',rgcP)
%    rgcVisualize('rgcTS',rgcP)
%    rgcVisualize('spikes',rgcP)
%    rgcVisualize('rf',rgcP)
%
%    rgcVisualize('spikesMeanMesh',rgcP,[3 1])
%    rgcVisualize('spikesMeanMesh',rgcP,[9 5])
%
%    rgcVisualize('cone voltages',rgcP)
%    rgcVisualize('cone voltages Mean',rgcP)
%    rgcVisualize('cone voltages Histogram',rgcP)
%
%    rgcVisualize('absorptions',rgcP)
%
% See also:  s_rgcHarmonic (illustrates calls)
%
% (c) Stanford Vista, Wandell, 2010

if notDefined('plotType'),error('Plot type undefined'); end
if notDefined('rgcP'),    error('rgcP required.'); end
if isempty(varargin), whichLayer = 1;
else                  whichLayer = varargin{1};
end

% number of simulated RGC layers
layer     = rgcP.get('layer',whichLayer);

% Lower case and remove spaces
plotType = stringFormat(plotType);

switch lower(plotType)

    % Input images
    case {'conevoltages','cvolts'}
        % rgcVisualize('cone voltages',rgcP);
        %
        % These are set by the cone adaptation model.
        % Someday this visualization should look like
        %   cones = rgcP.get('cones');
        %   conesVisualize(cones,'cone voltages');
        %
        cVolts = rgcP.get('cone voltages');
        % vSwing = pixelGet(sensorGet(rgcP.get('sensor'),'pixel'),'voltage swing');
        cVolts = 255*(cVolts/max(cVolts(:)));
        % Update call to mplayData ... see below.
        mplayData('cvolts',cVolts,layer);
    case {'cvoltsmean','conevoltagesmean' 'conevoltsmean'}
        cVolts = rgcP.get('cVolts');
        cVolts = mean(cVolts,3);
        cGrid  = rgcP.get('cone grid');
        imagesc(cGrid{1}, cGrid{2}, cVolts);         
        axis image
        title('Mean cone voltages')
        xlabel('cone position (microns)')
        ylabel('cone position (microns)')
        colormap gray
        colorbar
    
    case {'cvoltshistogram','conevoltageshistogram'}
        % The mean over time, and then a histogram across space
        cVolts = rgcP.get('cVolts');
        cVolts = mean(cVolts,3);
        hist(cVolts(:),50);

        % Absorptions
    case {'absorptions'}
        % rgcVisualize('absorptions',rgcP);
        absorptions = rgcP.get('absorptions');
        absorptions = 255*(absorptions - min(absorptions(:))/max(absorptions(:)));
        % Update call to mplayData ... see below.
        mp = mplay(absorptions);
        mFig = mp.hfig;
        set(mFig,'name','Cone absorptions')
        
    case {'absorptionsmean'}
        absorptions = rgcP.get('absorptions');
        absorptions = mean(absorptions,3);
        imagesc(absorptions);
    case {'absorptionshistogram'}
        % The mean over time, and then a histogram across space
        absorptions = rgcP.get('absorptions');
        absorptions = mean(absorptions,3);
        hist(absorptions(:),50);

        % Plots of the RGC properties
    case {'connectm'}
        error('NYI');
        % Do something better here.  This is just a print out, not a
        % visualization, of the connections.  And I don't understand the
        % arguments properly.  And we need an example.
        % connectm is a cell array, one cell for each RGC
        % Each cell contains two matrices.  The first matrix is a vector of
        % the other connected RGCs.
        % The second matrix are the time series (scaled by the connection
        % weight) of the coupling impulse response (cptr).
        % A visualization might be each point with a line connecting it to
        % its connected RGCs.
        if  isstruct(rgcData), connectM = rgcData.connectM;
        else connecMatrix = rgcData;
        end
        % if the option same layer is 1, then there could be 1 connecMatrix for
        % several layers.
        nMatrices = length(connectM);
        for im = 1:nMatrices
            middlePix = floor([imageSize(1) imageSize(2)]/2);
            middleIndex = sub2ind(imageSize,middlePix,middlePix);
            txt = strcat({'Layer '},int2str(im),{': '},'[',int2str(middlePix),{'] '},{'is connected to: '});
            connected = connectM{im}{middleIndex,1};
            for ii = 1:size(connected,1)
                index = connected(ii);
                [a b] = ind2sub(imageSize,index);
                txt = strcat(txt,{' ['}, int2str(a), {' '}, int2str(b), ']');
            end
            disp(txt{1});
        end

    case {'rf','rfimage','rfcenter','rfsurround','rfmesh','rfcentermesh','rfsurroundmesh'}
        % Plot the receptive field in this layer (combined center and
        % surround)
        % whichLayer = 1;
        % rgcVisualize('rf image',rgcP);
        % rgcVisualize('rf mesh',rgcP);
        % rgcVisualize('rf center',rgcP,whichLayer);

        
        % When we get the spatial RF, the get function accounts for the
        % scaling coefficients. This is NOT the case when we get the
        % temporal impulse response function. The return is just the
        % unscaled center and surround impulse response functions.
        vcNewGraphWin;
        if strfind(plotType,'center')
            rf = layer.get('RF center');
            tstr = 'RGC Receptive field center';
        elseif strfind(plotType,'surround')
            rf = layer.get('RF surround');
            tstr = 'RGC Receptive field surround';
        else
            rf = layer.get('RF sum');
            tstr = 'RGC Receptive field';
        end
        
        % This should be a get call.
        % Compute x,y axes in microns
        g = layerGet(layer,'rf grid','um');
        
        % See if you can figure out why Mathworks is annoying from this.
        if strfind(plotType,'mesh'),  mesh(g.X,g.Y,rf)
        else   imagesc(g.X(:),g.Y(:),rf), axis image; colorbar;
        end
        
        % Label
        xlabel('Position (um)')
        ylabel('Position (um)')
        title(tstr);

        % Attach to figure.
        udata.rf = rf;
        udata.x = g.X; udata.y = g.Y;
        set(gcf,'userdata',udata);

        % Temporal impulse response function
    case {'tirf','centertirf','surroundtirf','alltirf','intr'}
        % rgcVisualize('all tirf',rgcP)
        % rgcVisualize('center tirf',rgcP)
        % rgcVisualize('surround tirf',rgcP)
        % rgcVisualize(' tirf',rgcP)
        
        % When we get the spatial RF, the get function accounts for the
        % scaling coefficients. This is NOT the case when we get the
        % temporal impulse response function. The return is just the
        % unscaled center and surround impulse response functions.
        % The amplitude scaling is accounted for below.
        tIRF = layer.get('intr');
        dt = layer.get('dt');
        tSteps = dt*(1:size(tIRF,1));

        % Maybe there should a single call for getting scaled tIRF
        coefs = layer.get('rf coefs');
        tIRF  = tIRF*diag(coefs');
        
        vcNewGraphWin;
        % The TIRF is multipled by the sign of the spatial RF.
        % The spatial RF of the surround is opposte the sign of the center.
        % So, to see the TIRF here, we need to
        if strfind(plotType,'center')
            p = plot(tSteps,tIRF(:,1),'-k');            
            title(sprintf('Center Temporal Response Function: %s, [%1.2f %1.2f]', ...
                layer.get('name'),coefs(1),coefs(2)))
        elseif strfind(plotType,'surround')
            p = plot(tIRF(:,2),'-k');
            title(sprintf('Surround Temporal Response Function: %s, [%1.2f %1.2f]', ...
                layer.get('name'),coefs(1),coefs(2)))
        elseif strfind(plotType,'all')
            p = plot(tSteps,sum(tIRF,2),'-k', ...
                tSteps,tIRF(:,1),'g:', ...
                tSteps,tIRF(:,2),'r:');
            legend('Sum','Center','Surround')
            title(sprintf('Temporal Response Functions: %s, [%1.2f %1.2f]', ...
                layer.get('name'),coefs(1),coefs(2)))            
        else
            % Both tirf and intr plotted here
            p = plot(tSteps,tIRF(:,1) + tIRF(:,2),'-k');
            title(sprintf('Temporal Response Function (center + surround): %s, [%1.2f %1.2f]', ...
                layer.get('name'),coefs(1),coefs(2)))
            
        end

        % Adjust graph
        set(p,'linewidth',2);  xlabel('Time (ms)'); ylabel('Response (v)');

        
        grid on
        udata.tIRF = tIRF; udata.tSteps = tSteps;set(gcf,'userdata',udata);

    case {'feedbacktr','fbtr'}
        % rgcVisualize('feedback tr',rgcP)
        fbtr = layer.get('fbtr');
        t = rgcP.get('tr samples','ms');
        vcNewGraphWin;
        p = plot(t,fbtr);
        vThresh = layer.get('rgc volt thresh');
        hold on; 
        l = line([min(t),max(t)],[vThresh vThresh]);
        set(l,'LineWidth',2,'LineStyle','--','Color','r')
        hold off
        
        udata.fbtr = fbtr; udata.tSteps = t; udata.vThresh = vThresh;
        legend('Feedback','Voltage threshold')
        set(gcf,'userdata',udata);
        set(p,'linewidth',2);  xlabel('Time (ms)'); ylabel('Volts');
        title('Feedback'); grid on
        
    case {'couplingtr','cptr'}
        % rgcVisualize('coupling tr',rgcP)
        cptr = layer.get('cptr');
        t = rgcP.get('tr samples','ms');
        vcNewGraphWin; 
        p = plot(t,cptr);
        hold on;
        vThresh = layer.get('rgc volt thresh');
        l = line([min(t),max(t)],[vThresh vThresh]);
        set(l,'linewidth',2,'linestyle','--','color','r')
        hold off
        
        udata.fbtr = cptr; udata.tSteps = t;udata.vThresh = vThresh;
        legend('Coupling','Voltage threshold')
        set(gcf,'userdata',udata);
        set(p,'linewidth',2);  xlabel('Time (ms)'); ylabel('Volts');
        title('Coupling'); grid on
        
        % Plots of the RGC outputs
    case {'lints','lineartimeseries'}
        % rgcVisualize('Linear Time Series',rgcP,whichLayer);
        % rgcVisualize('Linear Time Series',rgcP);
        %
        % The linear time series response of each RGC is shown as a movie.

        linTS = layer.get('Linear TS');

        nT = size(linTS,2);
        gs = layer.get('grid size');
        linTS = reshape(linTS,[gs(1),gs(2),nT]);

        linTS = linTS - min(linTS(:));
        linTS = 255*linTS/max(linTS(:));
        % Update call to mplayData ... see below.
        mplayData('lints',linTS,layer);
        
    case {'componentts','componentimeseries'}
        % rgcVisualize('Component Time Series',rgcP,whichLayer);
        % rgcVisualize('Component Time Series',rgcP);
        %
        % The linear time series response as separate movies for each RGC
        % component (e.g., center and surround) shown as a movie

        componentTS = layer.get('component TS');
        nComponents = size(componentTS, 2);
        nT = size(componentTS{1},2);
        gs = layer.get('grid size');
        componentTSmovie = [];
        for cc = 1:nComponents
            tmp = reshape(componentTS{cc},[gs(1),gs(2),nT]);
            componentTSmovie = [componentTSmovie tmp];
        end
        componentTSmovie = componentTSmovie - min(componentTSmovie(:));
        componentTSmovie = 255*componentTSmovie/max(componentTSmovie(:));

        % Update call to mplayData ... see below.
        mplayData('componentts',componentTSmovie,layer)
        
    case {'lintsimage','lineartimeseriesimage' 'lintsmean'}
        % rgcVisualize('Linear Time Series Image',rgcP,whichLayer);
        % rgcVisualize('Linear Time Series Image',rgcP);
        %
        % The mean of the linear time series response of each RGC is shown as an image.

        linTS = layer.get('Linear TS');

        nT = size(linTS,2);
        gs = layer.get('grid size');
        linTS = reshape(linTS,[gs(1),gs(2),nT]);
        linTS = mean(double(linTS),3);
        
        vcNewGraphWin; imagesc(linTS);
        axis image equal off; colormap(gray);% colorbar;
       % title('Mean RGC linear time series (volts?)');


    case {'rgcts','finaltimeseries'}
        % rgcVisualize('RGC Time Series',rgcP,whichLayer);
        % rgcVisualize('RGC Time Series',rgcP);
        %
        % The time series response of each RGC is shown as a movie. Time
        % series is cone-driven + feedback + coupling

        rgcTS = layer.get('RGC TS');

        nT = size(rgcTS,2);
        gs = layer.get('grid size');
        rgcTS = reshape(rgcTS,[gs(1),gs(2),nT]);

        rgcTS = rgcTS - min(rgcTS(:));
        rgcTS = 255*rgcTS/max(rgcTS(:));
        
        % Update call to mplayData ... see below.
        mplayData('rgcts',rgcTS,layer);
       
    case {'rgctsmean','finaltimeseriesmean', 'rgctsimage'}
        % rgcVisualize('RGC Time Series Image',rgcP,whichLayer);
        % rgcVisualize('RGC Time Series Image',rgcP);
        %
        % The mean of the time series response of each RGC is shown as an
        % image. Responses are cone-driven + feedback + coupling

        rgcTS = layer.get('RGC TS');

        nT = size(rgcTS,2);
        gs = layer.get('grid size');
        rgcTS = reshape(rgcTS,[gs(1),gs(2),nT]);
        rgcTS = mean(double(rgcTS),3);
        vcNewGraphWin; imagesc(rgcTS);
        axis image equal off; colormap(gray); colorbar;
        title('Mean RGC time series (volts?)');


    case {'spikesmovie','spikes', 'spikests', 'spikets'}
        % The spikes as a time series.
        spikes = layer.get('spike time series');
        spikes = spikes*255;
        
        nT = size(spikes,2);
        gs = layer.get('grid size');
        spikes = reshape(spikes,[gs(1),gs(2),nT]);
        
        % Update call to mplayData ... see below.
        mplayData('spike time series',spikes,layer);
 
    case {'spikerate','spikerateimage', 'spikesmean', 'spiketsmean', 'spiketsimage', 'spikestsmean'}
        % rgcVisualize('Spike Rate Image',rgcP);
        % Average number of spikes over the entire stimulus
        % Could find out the stimulus duration and normalize per second.

        % GET THE STIMULUS DURATION !!!

        spikes = layer.get('spike time series');
        mSpikes = mean(double(spikes),3);
        vcNewGraphWin; imagesc(mSpikes);
        axis image equal off; colormap(gray); colorbar;
        title('Spike Rate (per frame)');

    case {'spikemesh','spikeratemesh'}
        % rgcVisualize('Spike Rate Mesh',rgcP);

        % Potentially blurred mesh of spike rate
        % Need to account for time properly here.
        spikes = layer.get('spike time series');
        mSpikes = mean(double(spikes),3);

        % Prepare to blur the mesh
        if isempty(varargin) || length(varargin) < 2, s = 7; w = 3;
        else s = varargin{2}(1); w = varargin{2}(2);
        end
        g = fspecial('gaussian', s, w);

        blurMSpikes = conv2(mSpikes,g,'same');
        vcNewGraphWin; mesh(blurMSpikes); colormap(cool)
        title('Spike Rate (per frame)');
    otherwise
        error('Unknown plotType %s\n',plotType);
end

return


%% Interface to mplay that the figure and controls the size
function mplayData(type,data,varargin)

switch stringFormat(type)
    case {'cvolts', 'conevoltages'}        
        mp = mplay(data); mFig = mp.hfig;
        set(mFig,'name', 'Cone Voltages')
        
    case {'lints','lineartimeseries'}
        if   isempty(varargin), error('layer needed'); 
        else layer = varargin{1};
        end
        
        mp = mplay(data); mFig = mp.hfig;
        rfCoefs = layer.get('rfcoeffs');       
        set(mFig,'name',sprintf('Linear TS: %s, [%1.2f %1.2f]', ...
            layer.get('name'),rfCoefs(1),rfCoefs(2)))
        
    case {'componentts'}
        if   isempty(varargin), error('layer needed');
        else layer = varargin{1};
        end

        rfCoefs = layer.get('rfcoeffs');
        mp = mplay(data); mFig = mp.hfig;
        set(mFig,'name',sprintf('Component TS: %s, [%1.2f %1.2f]', ...
            layer.get('name'),rfCoefs(1),rfCoefs(2)))

    case {'rgcts'}
        if isempty(varargin), error('layer needed');
        else layer = varargin{1};
        end
        
        mp = mplay(data); mFig = mp.hfig;
        rfCoefs = layer.get('rfcoeffs');
        set(mFig,'name',sprintf('RGC TS: %s, [%1.2f %1.2f]', ...
            layer.get('name'),rfCoefs(1),rfCoefs(2)))
        
    case {'spiketimeseries','spikets'}
        if isempty(varargin), error('layer needed');
        else layer = varargin{1};
        end

        mp = mplay(data); mFig = mp.hfig;
        rfCoefs = layer.get('rfcoeffs');
        set(mFig,'name',sprintf('Spike TS: %s, [%1.2f %1.2f]', ...
            layer.get('name'),rfCoefs(1),rfCoefs(2)))

    otherwise
        error('Unknown type %s\n',type);
end

% Make sure window size is reasonable.  If the data are small, the window
% is small.  So, we check and never let the window height be smaller than
% the width.
pos = get(mFig,'position');
if size(data,2) < pos(3)
    pos(4) = pos(3); set(mFig,'position',pos); 
end

return

