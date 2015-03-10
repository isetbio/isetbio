function res = layerGet(layer,fieldName,varargin)
% Implements the rgcLayer class get function 
%
%   res = layerGet(layer,fieldName,varargin)
%
% layer:     A layer stored in the class rgcParameters
% fieldName: string describing the field that we want
%
% This function is normally called by layer.get.  
%
% There are cases in which we use this function, called as
% layerGet(layer,param,varargin) directly.  We do this because we couldn't
% figure out how to make layer.get deal with the varargin that we use for
% setting units and other occasional varargins we will want.
%
%   res = layer.get(fieldName);
%
% The layer methods and parameters are defined in
%
%   methodsview('rgcLayer')
%   properties('rgcLayer')
%
% Parameter definitions for rgcP and rgcL are described in
% layerGet for the moment.  They will probably move elsewhere
% before long.
%
% Examples:
%   rgcP = rgcParameters;     % Initialize the parameters
%   rgcP.addOnLayer();        % Add a layer to the parameters
%   layer = rgcP.get('layer');
%   layer.get('type')
%   c = layer.get('center');
%
% Add aliases for the parameter in the function rgcMapParameterField.
%
% (c) Stanford Vista, 2010

if notDefined('layer'), error('layer object needed'); end
if notDefined('fieldName'), error('Field name needed'); end

% Spatial units, by default, are microns.
units = 'um';

% we map field name to a standard field
correctfield = rgcMapParameterField(fieldName);

% correctfield is a standard way of calling the field
switch(correctfield)
    case{'name'}
        % Book keeping
        res = layer.name;
        
    case {'type'}
        res = layer.type;
        
    case{'parent'}
        res = layer.parent;

        % Electrical parameters
    case {'vswing'}
        res = layer.vSwing;
        
        % Spatial parameters
    case{'cellspacing'}
        % Data are stored in um
        % layer.get('cell spacing','um')
        %
        if isempty(varargin), units = 'um';
        else units = varargin{1};
        end
        res = layer.cellSpacing;
        
        if ~strcmp(units,'um')
            res = (res/(10^6)) * ieUnitScaleFactor(units); 
        end

    case {'gridsize'}
        res = layer.gridSize;
        
    case {'cellgrid'}
        % These are the cell positions in microns with (0,0) at the center.
        % The positions are represented as two cell arrays, one for rows
        % and one for columns
        res = layer.cellGrid;
        
    case{'cellloc'}
        % These are the (x,y) positions of the cells in units of microns as
        % a matrix (nRow x nCol, 2).
        res = layer.cellLoc;
        
    case {'cellconepos'}
        % The cell positions with respect to the cone indices (integers)
        res = layer.cellConePos;
        
    case {'layercenter'}
        res = layer.layerCenter;
        
        % Temporal parameters
    case{'dt'}
        % Temporal step size in ms
        res = layer.dT;

    case{'trdur'}
        res = layer.trDur;
        
    case{'centertr'}
        % Should be a vector
        res = layer.centerTR;
        
    case{'surroundtr'}
        % Should be a vector
        res = layer.surroundTR;

%     case{'centerbdiv'}
%         res = layer.centerBDiv;
%         
%     case{'centerf'}
%         res = layer.centerF;
%         
%     case{'centernormf'}
%         res = layer.centerNormF;
%         
%     case{'centergamma'}
%         res = layer.centerGamma;
%         

%         
%     case{'surroundbdiv'}
%         res = layer.surroundBDiv;
%         
%     case{'surroundf'}
%         res = layer.surroundF;
%         
%     case{'surroundnormf'}
%         res = layer.surroundNormF;
%         
%     case{'surroundgamma'}
%         res = layer.surroundGamma;
%         
%     case{'gammaftr'}
%         res = layer.gammafTR;
%         
%     case{'gammabdiv'}
%         res = layer.gammaBDiv;
%         
%     case{'gammaf'}
%         res = layer.gammaF;
%         
%     case{'gammanormf'}
%         res = layer.gammaNormF;
%         
%     case{'gammagamma'}
%         res = layer.gammaGamma;
%         
%     case{'tshift'}
%         res = layer.tShift;
%         
%     case{'linf'}
%         res = layer.linF;
        
    case{'intr'}
        % Input (stimulus-driven) temporal response function
        res = layer.inTR;
        
        % Temporal and compute related
    case{'fbtr'}
        % Feedback temporal response
        res = layer.fbTR;
    
    case{'cptr'}
        % Coupling temporal response
        res = layer.cpTR;
        
        %% Spatial RF terms
    case{'rfcomponents','rf'}
        % Returns N-Dimensional matrix
        res = layer.RF;
    case{'rfsum','rfcomponentsum'}
        % Receptive field spatial structures of each component
        %
        % This got changed in June 2011 and edited again in Dec. 2011.  It
        % is confusing and potentially a problem.
        res = layer.RF;
        res = sum(res,3);
    case {'rfsupport'}
        % Spatial support of the whole receptive field.  
        % Default units: Microns.  
        % We want the support to be +/- 3SD for the rf. Hence the total
        % support should be 6 SD.
        
        if ~isempty(varargin), units = varargin{1}; end
        res = layerGet(layer,'rf std',units);
        res = 6*[max(res), max(res)];
        
    case {'rfgrid'}
        % layer.get('rf grid','um')
        % The (x,y) values of the RF support are returned as res.X and
        % res.Y.  The res.X and res.Y format is nice for
        % plotting.
        % The entire set of coordinates is [res.X(:), res.Y(:)];  
        
        if ~isempty(varargin), units = varargin{1}; end
        rfSupport = layerGet(layer,'rf support',units);
        
        rgcP = layer.get('parent');
        coneSpacing = rgcP.get('cone spacing',units);
        
        xPoints = 1:coneSpacing:rfSupport(2); xPoints = xPoints - mean(xPoints(:));
        yPoints = 1:coneSpacing:rfSupport(1); yPoints = yPoints - mean(yPoints(:));
        [X Y] = meshgrid(xPoints, yPoints);
        res.X = X;
        res.Y = Y;

    case{'rfcenter'}
        % Receptive field center spatially samples.
        % These are stored in units of ???
        res = layer.RF;
        res = res(:,:,1);
        
    case{'rfsurround'}
        % Receptive field surround spatially samples.
        % These are stored in units of ???
        res = layer.RF;
        res = res(:,:,2);
            case{'rfcov'}
        res = layer.RFcov;
        
    case {'rfstd'}
        % The units are 'um' by default
        % foo = layerGet(layer,'rf std',units);
        if ~isempty(varargin), units = varargin{1}; end
        
        cov = layer.RFcov;  % In units of um.
        lc = length(cov);   % Number of RF components, typically 2
        res = zeros(1,lc);
        for ii = 1:lc
            tmp = rgcExtractStdFromCovariance(cov{ii});
            res(ii) = tmp.std(1);
        end
        if ~strcmp(units,'um')
            res = (res/10^6) * ieUnitScaleFactor(units);
        end
        
    case {'rfellipse'}    
        % The default units appear to be microns. We could scale the std
        % units if we wanted.
        cov = layer.RFcov;  % In units of um.
        lc = length(cov);   % Number of RF components, typically 2
        res = cell(1,lc);
        for ii = 1:lc
            res{ii} = rgcExtractStdFromCovariance(cov{ii});
        end
        
    case{'rfscoeffs'}
        res = layer.RFsCoeffs;

    case {'coneweights'}
        res = layer.coneWeights;
        
    case {'cutoff'}
        res = layer.cutoff;
        
        % I think this might be some kind of parameter for the distance
        % function.
    case {'wrange'}
        res = layer.wRange;
        
        % Not sure what this is
    case {'nrange'}
        res = layer.nRange;
        
        % Voltage terms
    case {'rgcvoltthresh'}
        % This is a fraction of the RGC voltage swing 
        % Conceptually, this should be simplified.
        res = layer.rgcVoltThresh;
        
    case {'rgcvtimeseries'}
        % This is RGC time series after coupling and feedback, but before
        % thresholding for the spikes.
        % Name should be changed to rgcTS.
        res = layer.rgcvTimeSeries;
        
        % Spike time series properties
    case {'currentspkts','spkts'}
        % Spike time series (uint8) after RGC Time series is thresholded.
        % These are currently in the format of a (row,col,t) matrix. But
        % that might change.
        res = layer.currentSpkTS;
        
        % Linear time series properties
    case {'currentlints','lints'}
        % Linear time series after spatial and temporal convolution of the
        % cone absorptions by the RGC receptive field. This is prior to the
        % coupling and feedback processing in the RGC layer
        res = layer.currentLinTS;
    case {'maxlints'}
        res = max(layer.currentLinTS(:));
    case {'rangelints'}
        % layer.get('range lin ts')
        res = [min(layer.currentLinTS(:)),max(layer.currentLinTS(:))];
        
    case {'rfcomponentts'}
        % Linear signals after RF spatial and temporal filtering but before
        % combining across RF components (e.g., center and surround)
        % cell array 1 x number of RF components (typically 2) 
        res = layer.RFcomponentTS;
        
        % Connections
    case {'currentconnec'}
        res = layer.currentConnec;
        
        % Ugh.  Some calculation junk from GL
    case {'bordersize'}
        res = layer.getBorderSize();
        
    case {'overridesizedefault'}
        res = layer.overrideSizeDefault;
    
    case {'overridensize'}
        res = layer.overridenSize;
        
        % Probably shouldn't be here.
    case {'conespacing'}
        fprintf('cone spacing is an rgcP value')
        res = layer.parent.coneSpacing;

    otherwise
        error('Unknown field name');
end
