function [theFile, theRecord] = sdrPrebakedCircuitFile(mRGCMosaicFilename)
% Resolve a pre-baked ON-mRGC circuit filename to a local file
%
% Synopsis
%    [theFile, theRecord] = sdrPrebakedCircuitFile(mRGCMosaicFilename)
%
% Brief description
%   Pre-baked ON-midget-RGC circuits are published in the isetbio-mosaics
%   SDR deposit under normalized names. This function translates the
%   historical ISETBio filename that mRGCMosaic.loadPrebakedMosaic builds
%   into the deposited asset, downloading and verifying it when it is not
%   already cached.
%
%   A circuit present in the local ganglioncells/mosaics gallery wins, so a
%   locally synthesized mosaic is used without a network request.
%
%   Selection is by legacy filename rather than by the manifest's
%   eccentricity and size fields. Staging dropped the sign of the x
%   eccentricity, so the published metadata cannot distinguish a circuit at
%   x = -2 degrees from one at x = +2 degrees. See sdrLegacyFilenameMap.
%
% Inputs
%   mRGCMosaicFilename - Historical ISETBio circuit filename
%
% Returns
%   theFile   - Full path to a local, verified circuit MAT-file, or empty
%               when the deposit has no matching asset
%   theRecord - The legacy-map record used, or empty
%
% See also
%   mRGCMosaic.loadPrebakedMosaic, mRGCMosaic.listPrebakedMosaics,
%   sdrLegacyFilenameMap, sdrMosaicRecords, sdrMosaicFetch

theFile = '';
theRecord = [];

% A locally synthesized circuit takes precedence over the deposited copy.
[~, prebakedMRGCMosaicDir] = mRGCMosaic.listPrebakedMosaics();
localFile = fullfile(prebakedMRGCMosaicDir, mRGCMosaicFilename);
if isfile(localFile)
    theFile = localFile;
    return;
end

mapRecords = sdrLegacyFilenameMap('mrgc_on_circuits');
mapIndex = find(strcmp({mapRecords.legacy_filename}, mRGCMosaicFilename), 1);
if isempty(mapIndex)
    return;
end
theRecord = mapRecords(mapIndex);

% Select the manifest record by path, then let sdrMosaicFetch verify it.
sdrRecords = sdrMosaicRecords('mrgc_on_circuits');
sdrIndex = find(strcmp({sdrRecords.path}, theRecord.sdr_path), 1);
if isempty(sdrIndex)
    error('sdrPrebakedCircuitFile:Manifest', ...
        ['The legacy filename map names %s, but the published deposit ', ...
        'has no such asset.'], theRecord.sdr_path);
end

theFile = sdrMosaicFetch('mrgc_on_circuits', sdrRecords(sdrIndex));

end
