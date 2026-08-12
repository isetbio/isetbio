function localFile = sdrMosaicFetch(collectionName, record)
% Return a verified local copy of one isetbio-mosaics asset
%
% Synopsis
%    localFile = sdrMosaicFetch(collectionName, record)
%
% Brief description
%   Resolve one manifest record to a file in the local cache, downloading it
%   from the published SDR deposit when necessary. The asset is written to a
%   temporary sibling, checked against the manifest byte count and SHA-256,
%   and only then renamed into its mirrored cache path. A cached file is
%   re-verified on every call, so a corrupt cache entry is replaced rather
%   than returned.
%
% Inputs
%   collectionName - Collection holding the asset, such as 'cone_mosaics'
%   record         - One manifest record from sdrMosaicRecords
%
% Returns
%   localFile - Full path to the verified local copy
%
% See also
%   sdrMosaicCacheRoot, sdrMosaicRecords, mosaicLoad, ieWebGet

cacheRoot = sdrMosaicCacheRoot();
pathParts = split(record.path, '/');
localFile = fullfile(cacheRoot, collectionName, pathParts{:});

if isfile(localFile) && sdrMosaicVerify(localFile, record)
    return;
end

% Either there is no cached copy or the cached copy no longer matches the
% manifest. Remove the bad file so a failed download cannot leave a stale
% asset in place.
if isfile(localFile), delete(localFile); end

temporaryRoot = tempname(cacheRoot);
mkdir(temporaryRoot);
cleanup = onCleanup(@() localRemoveFolder(temporaryRoot));

temporaryFile = ieWebGet('depositname', 'isetbio-mosaics', ...
    'depositfile', [collectionName '/' record.path], ...
    'downloaddir', temporaryRoot, 'unzip', false);

if ~sdrMosaicVerify(temporaryFile, record)
    error('sdrMosaicFetch:Checksum', ...
        ['Downloaded %s asset ''%s'' failed its manifest checksum or ', ...
        'byte-count check.'], collectionName, record.path);
end

localDir = fileparts(localFile);
if ~isfolder(localDir), mkdir(localDir); end
movefile(temporaryFile, localFile, 'f');

clear cleanup

end

function localRemoveFolder(folder)
if isfolder(folder), rmdir(folder, 's'); end
end
