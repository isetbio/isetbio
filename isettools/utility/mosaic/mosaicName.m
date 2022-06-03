function fname = mosaicName(sizeDegs,positionDegs)
% Return a standard mosaic name for data/cones file
%
% Mosaics with these names are stored in a small library in
% isettools/data/cones
%
% See also
%   mosaicLoad

fname = sprintf('cmosaic_%.1f-%.1f_%.1f-%.1f.mat',sizeDegs,positionDegs);
fname = fullfile(isetRootPath,'data','cones',fname);

end