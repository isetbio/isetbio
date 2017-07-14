function hdl = plot(obj, pType, varargin)
% Plot data from a bipolar mosaic object
% 
%    hdl = bp.plot(plotType, varargin)
%
% Plot types
%   response center
%   response surround
%   response image
%   response movie
%   response time series
%   spatial rf
%   mosaic
%
% Optional parameters-value pairs
%   gamma - controls image display
%   pos   - position
%
% Examples:
%   (run s_LayersTest).
%   bpMosaics = bpL.mosaic;
%   bpMosaics{1}.plot('spatial rf')
%   bpMosaics{1}.plot('mosaic');
%   bpMosaics{1}.plot('response center');
%   bpMosaics{1}.plot('response time series','pos',[5 5]);
%   bpMosaics{1}.plot('response image');
%   bpMosaics{1}.plot('response movie');
%
% 5/2016 JRG,BW (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName  = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'responsetimeseries','responsecenter','responsesurround',...
    'movieresponse','responseimage', ...
    'spatialrf','mosaic','responsemovie'};
p.addRequired('pType',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

p.addParameter('gamma',1,@isscalar);
p.addParameter('pos',[],@isvector);

% Parse pType
p.parse(pType,varargin{:}); 
% Parameters are pulled out in the case statements, below.

%% Create window
hdl = gcf;%vcNewGraphWin([],'upperLeftBig');

sz = size(obj.responseCenter);

% Programming:
% We need to get the units of time from the object, not as per below.

% Options
switch ieParamFormat(pType)
    case 'spatialrf'
        % bp.plot('spatial rf')
        % Not a good idea yet, because the rf is usually an impulse these
        % days.  May change in the future.
        srf = obj.sRFcenter - obj.sRFsurround;
        sz = size(srf); 
        if isequal(sz,[1, 1])
            disp('spatial rf is an impulse');
            return;
        end
        
        x = (1:sz(2)) - mean(1:sz(2));    
        y = (1:sz(1)) - mean(1:sz(1)); 
        surf(x,y,srf); colormap(parula);
        xlabel('Cone position re: center'); zlabel('Responsivity')
        
    case {'mosaic'}
        % bp.plot('mosaic') - Shows RF array
        % Get contour lines for mosaic RFs
                
        % Oddly, the center is (row,col)
        center = obj.cellLocation;  % um w.r.t. center of image
        center = reshape(center,[size(obj.cellLocation,1)*size(obj.cellLocation,2),2]);
        radius = 1e6*.5*obj.patchSize/(size(obj.cellLocation,1));
        ellipseMatrix = [1 1 0];        
        ieShape('ellipse','center',center,...
            'radius',0.5*radius,...
            'ellipseParameters',ellipseMatrix,...
            'color','b');
        
        % Sets the axis limits
        set(gca,...
            'xlim',[min(center(:,2)) - 3*radius, max(center(:,2)) + 3*radius],...
            'ylim',[min(center(:,1)) - 3*radius, max(center(:,1)) + 3*radius]);
        xlabel(sprintf('Distance (\\mum)'),'fontsize',14);
        ylabel(sprintf('Distance (\\mum)'),'fontsize',14);
            
    case{'responsecenter'}
        % bp.plot('response center')
        responseRS = reshape(obj.responseCenter,sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'responsesurround'}
        % bp.plot('response surround')
        responseRS = reshape(obj.responseSurround,sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'responsetimeseries'}
        % bp.plot('response time series','pos',...)
        %
        % Open a new window and show the time series, but not in any
        % particular organized way. Only up to 200 samples are shown.
        % Random draws.  We should allow this to be controlled.
        %
        pos = p.Results.pos;
        if isempty(pos)
            nCells = sz(1)*sz(2);
            response = reshape(obj.get('response'),nCells,sz(3));    
        else
            response = squeeze(obj.responseCenter(pos(1),pos(2),:) - obj.responseSurround(pos(1),pos(2),:));
        end
        tSamples = obj.timeStep*(1:size(obj.responseCenter,3));
        
        vcNewGraphWin;
        if isempty(pos) && nCells >= 100
            % 100 randomly sampled positions
            cSamples = randi(nCells,[100,1]);
            plot(tSamples,response(cSamples,:));
            title('Bipolar current (100 samples)');
        elseif isempty(pos)
            % All of the spatial samples
            plot(tSamples,response);
            title('Bipolar current');
        else
            % From a single point
            plot(tSamples,response)
            title(sprintf('Position %d,%d',pos(1),pos(2)));
        end
        
        xlabel(sprintf('Time (sec, \\Deltat %.0f ms)',obj.timeStep*1e3));
        ylabel('Current (AU)');
        grid on
        
    case{'responseimage'}
        % bp.plot('response')
        response = (obj.get('response'));
        patchSizeUM = obj.patchSize*1e6;
        dx = patchSizeUM/size(response,1);  % Step in microns
        samples = 1:dx:size(response,1);
        samples = samples - mean(samples(:));
        if p.Results.gamma ~= 1, img = mean(response,3).^p.Results.gamma;
        else                     img = mean(response,3);
        end
        
        imagesc(samples,samples,img);
        axis image; colormap(gray(256)); title('Current (a.u.)');
        xlabel(sprintf('Cell position (\\mum)'));
        
    case{'responsemovie'}
        % Pass the varargin along
        if ~isempty(varargin) && length(varargin) == 1
            % Params are coded in a single struct
            varargin{1}.hf = hdl;
            if p.Results.gamma ~= 1
                varargin{1}.gamma = p.Results.gamma;
            end
            ieMovie(obj.get('response'),varargin{:});
        else
            % Params are coded in param/value
            varargin{end+1} = 'gamma'; varargin{end+1} = p.Results.gamma;
            ieMovie(obj.get('response'),'hf',hdl,varargin{:});
        end
end

end
