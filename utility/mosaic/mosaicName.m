function fname = mosaicName(sizeDegs,positionDegs)
% Return a standard mosaic name for data/cones file
%
%   cmosaic_sizedegs_positiondegs
%
% For example, 1 deg on a size at the central fovea
%
%   cmosaic_1.0-1.0_0.0-0.0.mat
%
% The cmosaics are stored in a small library in isettools/data/cones
% (isetbioDataPath)
%
% See also
%   mosaicLoad

fname = sprintf('cmosaic_%.1f-%.1f_%.1f-%.1f.mat',sizeDegs,positionDegs);
fname = fullfile(isetbioDataPath,'cones',fname);

end