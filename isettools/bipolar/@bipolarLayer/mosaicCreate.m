function mosaicCreate(obj,cellType, varargin)
% MOSAICCREATE - Bipolar layer method to create a mosaic of bipolar cells
%
%   bpL = @bipolarLayer.mosaicCreate('cellType',varargin,...);
%
% Wandell, ISETBIO Team, 2017

%% Set the sell type and pass in the layer object which will become parent

p = inputParser;
p.KeepUnmatched = true;

% Check the validity in the mosaic class function.  We should make the
% valid mosaics part of the 'layer' object and check here.
p.addRequired('cellType',@ischar);
p.addParameter('nMosaic',[],@isinteger);

p.parse(cellType,varargin{:});

% Set up the cellType and parent parameters
params.cellType = p.Results.cellType;
params.parent   = obj;

% Create the mosaic as the next one in the list.  Maybe we should have a
% parameter to set this to a specific number.
nMosaic = p.Results.nMosaic;

if isempty(nMosaic)
    % Add it to the end
    obj.mosaic{end+1} = bipolarMosaic(obj.input, params);
else
    % Use the specified slot
    obj.mosaic{nMosaic} = bipolarMosaic(obj.input, params);
end

end
