function records = sdrMosaicRecords(collectionName)
% Return the manifest records for one isetbio-mosaics collection
%
% Synopsis
%    records = sdrMosaicRecords(collectionName)
%
% Brief description
%   Read the published manifest for a collection of the isetbio-mosaics SDR
%   deposit, downloading the top-level and collection manifests into the
%   local cache the first time they are needed.
%
%   Every asset in the deposit is described by one record. A loader selects
%   a record by its metadata and passes it to sdrMosaicFetch; loaders must
%   not reconstruct a remote filename from MATLAB arguments.
%
% Inputs
%   collectionName - One of 'cone_lattices', 'cone_mosaics',
%                    'mrgc_lattices', or 'mrgc_on_circuits'
%
% Returns
%   records - Struct array of manifest records. Every record has path,
%             bytes, sha256, kind, and eye fields; the remaining fields are
%             collection specific.
%
% See also
%   sdrMosaicCacheRoot, sdrMosaicFetch, mosaicLoad, ieWebGet

cacheRoot = sdrMosaicCacheRoot();

topManifest = localFetchManifest(cacheRoot, 'manifest.json');
topData = jsondecode(fileread(topManifest));

collectionIndex = find(strcmp({topData.collections.name}, collectionName), 1);
if isempty(collectionIndex)
    error('sdrMosaicRecords:UnknownCollection', ...
        'The SDR manifest has no %s collection. Valid collections are: %s.', ...
        collectionName, strjoin({topData.collections.name}, ', '));
end
theCollection = topData.collections(collectionIndex);

collectionManifest = localFetchManifest(cacheRoot, theCollection.manifest);
collectionData = jsondecode(fileread(collectionManifest));
records = collectionData.files;

% The top-level manifest states how many assets each collection holds. A
% count mismatch means a truncated or stale cached manifest, so replace it
% once before giving up.
if numel(records) ~= theCollection.file_count
    delete(collectionManifest);
    collectionManifest = localFetchManifest(cacheRoot, theCollection.manifest);
    collectionData = jsondecode(fileread(collectionManifest));
    records = collectionData.files;
    if numel(records) ~= theCollection.file_count
        error('sdrMosaicRecords:Manifest', ...
            'The %s manifest lists %d files but the deposit declares %d.', ...
            collectionName, numel(records), theCollection.file_count);
    end
end

end

function localFile = localFetchManifest(cacheRoot, relativePath)
% Download one manifest into the mirrored cache if it is not already there.
pathParts = split(relativePath, '/');
localFile = fullfile(cacheRoot, pathParts{:});
if ~isfile(localFile)
    ieWebGet('depositname', 'isetbio-mosaics', 'depositfile', relativePath, ...
        'downloaddir', cacheRoot, 'unzip', false);
end
end
