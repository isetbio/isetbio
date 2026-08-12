function theMap = sdrLegacyFilenameMap(collectionName)
% Legacy ISETBio filenames for the isetbio-mosaics SDR assets
%
% Synopsis
%    theMap  = sdrLegacyFilenameMap()
%    records = sdrLegacyFilenameMap(collectionName)
%
% Brief description
%   Return the map between the filenames ISETBio used before the SDR
%   migration and the paths of the deposited assets. Every pairing was
%   established by comparing the SHA-256 of a repository copy with the
%   published manifest, so the map is evidence rather than a naming
%   convention.
%
%   The map exists because the published manifests do not yet carry a
%   legacy_filename field, and because the mrgc_on_circuits manifest lost
%   the sign of the x eccentricity during staging: two circuits at x = -2
%   and x = +2 degrees both appear as ecc-x2-y0 with eccentricity_degrees
%   [2 0]. A circuit therefore cannot be selected from the published
%   metadata alone, and the legacy filename is used as the lookup key.
%
%   Five deposited circuits never had a repository copy. Their records
%   carry an empty legacy_filename plus a candidate_legacy_filename and a
%   candidate_confidence, and must not be treated as verified.
%
% Inputs
%   collectionName - Optional. One of 'cone_lattices', 'cone_mosaics',
%                    'mrgc_lattices', or 'mrgc_on_circuits'. When given,
%                    only that collection's records are returned.
%
% Returns
%   theMap - The whole decoded map, or one collection's record array
%
% See also
%   sdrMosaicRecords, sdrMosaicFetch, mosaicLoad

persistent cachedMap

if isempty(cachedMap)
    mapFile = fullfile(isetbioDataPath, 'sdr', ...
        'isetbio-mosaics-legacy-filenames.json');
    if ~isfile(mapFile)
        error('sdrLegacyFilenameMap:Missing', ...
            'The legacy filename map is missing from %s.', mapFile);
    end
    cachedMap = jsondecode(fileread(mapFile));
end

theMap = cachedMap;

if nargin > 0
    if ~isfield(theMap.collections, collectionName)
        error('sdrLegacyFilenameMap:UnknownCollection', ...
            'The legacy filename map has no %s collection.', collectionName);
    end
    theMap = theMap.collections.(collectionName).records;
end

end
