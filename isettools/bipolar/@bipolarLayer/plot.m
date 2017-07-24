function hdl = plot(obj, pType, varargin)
% Deprecate - plots should be called based on the mosaic, I think.
%
% Plot the values from one of the mosaics of a bipolarLayer object
% 
%    hdl = bipolarLayer.plot(plotType,...)
%
% Required Inputs
%  pType -  The plot type.  Use @bipolarLayer.plot('help','nMosaic',1) to
%           see implemented plot types.
%  
% Parameter-Key Inputs
%  nMosaic  - Mosaic number to plot
%  gamma    
%  pos
%  newWindow
%
% For most plot types this routine calls rgcMosaic.plot. The selected
% mosaic from the layer is based on the integer nMosaic.
%
% We plan to implement cases for the layer that compare across several
% mosaics.
%
% Examples:
%   bipolarLayer.plot('mosaic','nMosaic',1);
%   bipolarLayer.plot('movie response','nMosaic',1);
%
% JRG/BW (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName  = mfilename;
p.KeepUnmatched = true;

p.addRequired('pType',@ischar);

p.addParameter('gamma',1,@isscalar);
p.addParameter('pos',[1,1],@isvector);    % If not passed, then assume 1,1
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
    disp('No layer plots implemented yet.')
    return;
else
    % A plot based on one mosaic.  Call the bipolarMosaic.plot funciton.
    hdl = obj.mosaic{nMosaic}.plot(pType,p.Results);
end

end


