function hdl = plot(obj, pType, varargin)
% Plot the values from the bipolar object
% 
%    hdl = bp.plot(parameter)
%
% Plot types
%   response center
%   response surround
%   response
%   movie response 
%   Time series
% Optional parameters
%   gamma for image display
%
% Examples:
%
% 
% 5/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName  = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'responsetimeseries','responseCenter','responseSurround',...
    'movieresponse','responseimage', ...
    'spatialrf','mosaic'};
p.addRequired('pType',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

p.addParameter('gamma',1,@isscalar);
p.addParameter('pos',[],@isvector);

% Parse pType
p.parse(pType,varargin{:}); 



%% Create window
hdl = gcf;%vcNewGraphWin([],'upperLeftBig');

sz = size(obj.responseCenter);

% Programming:
% We need to get the units of time from the object, not as per below.

% Options
switch ieParamFormat(pType)
    case 'spatialrf'
        % bp.plot('spatial rf')
        srf = obj.sRFcenter - obj.sRFsurround;
        sz = size(srf); 
        x = (1:sz(2)) - mean(1:sz(2));    
        y = (1:sz(1)) - mean(1:sz(1)); 
        surf(x,y,srf); colormap(parula);
        xlabel('Cone position re: center'); zlabel('Responsivity')
    case {'mosaic'}
        % bp.plot('mosaic') - Shows RF array
        % At some point we will allow the sRF parameters to vary across the
        % array.
        % See irPlot.m for an example.
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
        
        %
        %         spatialRFcontours = plotContours(obj.mosaic{cellTypeInd}, obj.size, obj.col);
        %         figure(hf);
        %         % Subplot if more than one mosaic is being plotted
        %         if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
        %
        %         nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
        %
        %         % Convert RGC position to distance
        %         patchSizeX = obj.size;
        %         numberCellsX = obj.col;
        %         umPerCell = 1e0*patchSizeX/numberCellsX;
        %
        %         cmap = parula(16);
        %         for xcell = 1:nCells(1)
        %             for ycell = 1:nCells(2)
        %                 hold on;
        %                 % center
        %                 plot(umPerCell*spatialRFcontours{xcell,ycell,1}(1,2:end),...
        %                     umPerCell*spatialRFcontours{xcell,ycell,1}(2,2:end),...
        %                     'color','r','linewidth',3);%cmap(cellTypeInd,:));
        %                 hold on;
        %                 % surround
        %                 plot(umPerCell*spatialRFcontours{xcell,ycell,2}(1,2:end),...
        %                     umPerCell*spatialRFcontours{xcell,ycell,2}(2,2:end),...
        %                     'color','m','linewidth',3);%cmap(cellTypeInd+8,:));
        %             end
        %         end
        %         axis equal
        %         title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
        %         xlabel(sprintf('Distance (m)'),'fontsize',14);
        %         ylabel(sprintf('Distance (m)'),'fontsize',14);
            
            
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
        axis image; colormap(gray(256));
        xlabel(sprintf('Cell position (\\mum)'));
        
    case{'responsemovie','movieresponse'}
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
