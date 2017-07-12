function hdl = plot(obj, pType, varargin)
% Plot the values from one of the mosaics of a bipolarLayer object
% 
%    hdl = bipolarLayer.plot(plotType,...)
%
% Required Inputs
%  plotType - Plot type {response center, response surround, response movie,
%             response time series}
%  
% Parameter-Key Inputs
%  nMosaic  - Mosaic to use for plotting
%  gamma
%  pos
%  newWindow
%
% For many parameters, this routine calls rgcMosaic.plot. In that case the
% selected mosaic is based on the integer nMosaic.
%
% We will develop some cases in which we compare across mosaics. In that
% case nMosaic will be 0.
%
% Optional parameters
%   gamma for image display
%
% Examples:
%   bipolarLayer.plot('mosaic','nMosaic',1);
%   bipolarLayer.plot('movie response','nMosaic',1);
%
% 5/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName  = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.  This list is pretty small.  Hard to
% maintain. (BW).
allowableFields = {...
    'responsetimeseries','responseCenter','responseSurround',...
    'movieresponse','responseimage', ...
    'spatialrf','mosaic'};
p.addRequired('pType',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

p.addParameter('gamma',1,@isscalar);
p.addParameter('pos',[],@isvector);
p.addParameter('newWindow',false,@islogical);

% Will this be one of the mosaics, or use multiple mosaics (0, or maybe a
% vector?)
p.addParameter('nMosaic',0,@(x)(isscalar(x) && (x >= 0)));

% Parse pType
p.parse(pType,varargin{:}); 

nMosaic = p.Results.nMosaic;
if nMosaic > length(obj.mosaic)
    error('nMosaic (%d) exceeds number of bipolar layers (%d)\n',...
        nMosaic,length(obj.mosaic));
end

%% Account for parameters
if p.Results.newWindow; hdl = vcNewGraphWin; end

if nMosaic == 0
    % A plot that uses more than one mosaic.  What the layer is for.
    % We need to make some stuff up for here.
    disp('Whole layer plot options are NYI')
    return;
else
    % A plot based on one mosaic.  Call the bipolarMosaic.plot funciton.
    hdl = obj.mosaic{nMosaic}.plot(pType);
end

end


