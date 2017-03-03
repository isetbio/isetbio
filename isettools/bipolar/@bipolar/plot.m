function hdl = plot(obj, pType,varargin)
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
% 
% Examples:
%
% 
% 5/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'response','responseCenter','responseSurround',...
    'movieresponse', ...
    'spatialrf','mosaic'};
p.addRequired('pType',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

% Parse pType
p.parse(pType,varargin{:}); 

%% Create window
hdl = vcNewGraphWin([],'upperLeftBig');

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
        
        spatialRFcontours = plotContours(obj.mosaic{cellTypeInd}, obj.size, obj.col);
        figure(hf);
        % Subplot if more than one mosaic is being plotted
        if length(mosaicIndices)>1; subplot(ceil(length(mosaicIndices)/2),2,cellTypeInd); end;
        
        nCells = obj.mosaic{cellTypeInd}.get('mosaic size');
        
        % Convert RGC position to distance
        patchSizeX = obj.size;
        numberCellsX = obj.col;
        umPerCell = 1e0*patchSizeX/numberCellsX;
        
        cmap = parula(16);
        for xcell = 1:nCells(1)
            for ycell = 1:nCells(2)
                hold on;
                % center
                plot(umPerCell*spatialRFcontours{xcell,ycell,1}(1,2:end),...
                    umPerCell*spatialRFcontours{xcell,ycell,1}(2,2:end),...
                    'color','r','linewidth',3);%cmap(cellTypeInd,:));
                hold on;
                % surround
                plot(umPerCell*spatialRFcontours{xcell,ycell,2}(1,2:end),...
                    umPerCell*spatialRFcontours{xcell,ycell,2}(2,2:end),...
                    'color','m','linewidth',3);%cmap(cellTypeInd+8,:));
            end
        end
        axis equal
        title(sprintf('%s',obj.mosaic{cellTypeInd}.cellType),'fontsize',14);
        xlabel(sprintf('Distance (m)'),'fontsize',14);
        ylabel(sprintf('Distance (m)'),'fontsize',14);
            
            
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
        
    case{'response'}
        % bp.plot('response')
        response = reshape(obj.get('response'),sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),response);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'movieresponse'}
        % Pass the varargin along
        if ~isempty(varargin) && length(varargin) == 1
            % Params are coded in a single struct
            varargin{1}.hf = hdl;
            ieMovie(obj.get('response'),varargin{:});
        else
            % List of params
            r = obj.get('response');
            ieMovie(r,'hf',hdl,varargin{:});
        end
end

% set(gca,'fontsize',16);

end
