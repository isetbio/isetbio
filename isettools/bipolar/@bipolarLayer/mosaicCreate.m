function mosaicCreate(obj,cellType, varargin)
% MOSAICCREATE - Create a mosaic of bipolar cells
%
%   bpL = @bipolarLayer.mosaicCreate('cellType',...);
%
%
% Wandell, ISETBIO Team, 2017

%%
p = inputParser;
p.KeepUnmatched = true;

% We check the validity in the mosaic class function 
addRequired(p, 'cellType',@ischar);

p.parse(cellType,varargin{:});

%%
% bpMosaicParams.cellType = p.Results.cellType;
obj.mosaic{end+1} = bipolarMosaic(obj.input, p.Results);

end
