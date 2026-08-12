function cacheRoot = sdrMosaicCacheRoot()
% Root of the local cache mirroring the isetbio-mosaics SDR deposit
%
% Synopsis
%    cacheRoot = sdrMosaicCacheRoot()
%
% Brief description
%   Assets downloaded from the published isetbio-mosaics deposit are cached
%   here, preserving the deposited directory hierarchy:
%
%      <cacheRoot>/<collection>/<eye>/<file>.mat
%
%   The directory is excluded from Git by data/sdr/.gitignore, so a cached
%   asset is never committed. Delete the directory to force a re-download.
%
% Inputs
%   None
%
% Returns
%   cacheRoot - Full path to the local isetbio-mosaics cache root
%
% See also
%   sdrMosaicRecords, sdrMosaicFetch, mosaicLoad, ieWebGet

cacheRoot = fullfile(isetbioRootPath, 'data', 'sdr', 'isetbio-mosaics');

end
