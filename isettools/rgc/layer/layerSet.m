function layer = layerSet(layer,fieldName,val)
% implements the set function for the rgcLayer class
%
% layer = layerSet(layer,fieldName,val)
%
% layer:     the parameter object (from the class rgcParameters)
% fieldName: string describing the field that we want
% val:     val to be set to
%
% This function is equivalent to
%       layer.set(fieldName, val);
% which is probably the way it should be called
%
% Example:
%   If you feel like an alias is missing for a field, feel free to add it
%   in rgcMapParameterField. But this function only uses the standard name. 
%
% Stanford, 2011

if notDefined('layer'),  error('layer object needed'); end
if notDefined('fieldName'), error('Field name needed'); end

% we map field name to a standard field
correctfield = rgcMapParameterField(stringFormat(fieldName));

% correctfield is a standard way of calling the field
switch(correctfield)
    case{'name'}
        layer.name = val;
        
    case{'cellspacing'}
        layer.cellSpacing = val;
        
    case {'gridsize'}
        layer.gridSize = val;
        
    case {'cellgrid'}
        layer.cellGrid = val;
        
    case{'cellloc'}
        layer.cellLoc = val;
        
    case{'parent'}
        layer.parent = val;
        
    case{'rfcov'}
        % This is a cell array of two covariance matrices
        % {[varC 0; 0 varC] [varS 0; 0 varS]}
        % These are the standard deviation of the center and the surround,
        % each of which can be oriented and/or anisotropic.  
        layer.RFcov = val;
        
    case {'rfstd'}
        % Convert a vector of standard deviations and orientations into a
        % cell array of covariance matrices.
        %   val{1} = {1,  3,  pi/2};
        %   layer.set('rfstd',val);
        %
        lc = length(val);  % How many?
        rfCov = cell(1,lc);  % Create the cell array of covariance matrices
        
        for ii = 1:lc
            rfCov{ii} = rgcCovarianceMatrix((val{ii}(1)),(val{ii}(2)),val{ii}(3));
        end
        % Set the cell array
        layer.RFcov = rfCov;
        
       
    case{'rfscoeffs'}
        layer.RFsCoeffs = val;
    
    case {'coneweights'}       
        layer.coneWeights = val;
        
    case{'trdur'}
        layer.trDur = val;
        
    case{'centertr'}
        layer.centerTR = val;
        
    case{'centerbdiv'}
        layer.centerBDiv = val;
        
    case{'centerf'}
        % 
        layer.centerF = val;
        
    case{'centernormf'}
        layer.centerNormF = val;
        
    case{'centergamma'}
        layer.centerGamma = val;
        
    case{'surroundtr'}
        layer.surroundTR = val;
        
    case{'surroundbdiv'}
        layer.surroundBDiv = val;
        
    case{'surroundf'}
        layer.surroundF = val;
        
    case{'surroundnormf'}
        layer.surroundNormF = val;
        
    case{'surroundgamma'}
        layer.surroundGamma = val;
        
    case{'gammaftr'}
        layer.gammafTR = val;
        
    case{'gammabdiv'}
        layer.gammaBDiv = val;
        
    case{'gammaf'}
        layer.gammaF = val;
        
    case{'gammanormf'}
        layer.gammaNormF = val;
        
    case{'gammagamma'}
        layer.gammaGamma = val;
        
    case{'dt'}
        layer.dT = val;
        
    case{'intr'}
        layer.inTR = val;
        
    case{'fbtr'}
        layer.fbTR = val;
        
    case{'cptr'}
        layer.cpTR = val;
        
        % Electrical properties
    case {'vswing'}
        layer.vSwing = val;
        
    case {'rgcvoltthresh'}
        layer.rgcVoltThresh = val;
        
           % Network properties
    case {'cutoff'}
        layer.cutoff = val;
        
    case {'wrange'}
        layer.wRange = val;
        
    case {'nrange'}
        layer.nRange = val; 
        % Spatial properties
        
    case{'rf'}
        % What is this?
        layer.RF = val;
        
    case {'layercenter'}
        layer.layerCenter = val;
        
    case {'rgcvtimeseries'}
        % The RGC time series
        layer.rgcvTimeSeries = val;
        
    case {'currentspkts'}
        % The RGC spike series
        layer.currentSpkTS = val;
        
    case {'currentlints'}
        % Linear signals after RF spatial and temporal filtering
        layer.currentLinTS = val;
        
    case {'rfcomponentts'}
        % Linear signals after RF spatial and temporal filtering but before
        % combining across RF components (e.g., center and surround)
        % cell array 1 x number of RF components (typically 2) 
        layer.RFcomponentTS = val;
       
    case {'currentconnec'}
        layer.currentConnec = val;
        
    case {'hascoupling'}
        layer.hasCoupling = val;        
        
    case {'hasfeedback'}
        layer.hasFeedback = val;
        
    otherwise
        error('Unknown field name or unauthorized operation');
end
