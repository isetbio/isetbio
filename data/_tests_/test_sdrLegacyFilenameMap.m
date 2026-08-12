function tests = test_sdrLegacyFilenameMap()
% Offline checks on the isetbio-mosaics legacy-filename map
%
% The map records the pre-SDR ISETBio filename for every deposited asset.
% It is the lookup table used while the published manifests still lack a
% legacy_filename field, so its integrity is checked without a network.
%
% The parity check below is the guard used before deleting a bundled asset:
% a repository copy must be byte-identical to the deposited copy the map
% names. Assets already migrated to the SDR have no repository copy left, so
% those records are skipped rather than failed.
%
% See also
%   sdrLegacyFilenameMap, sdrMosaicRecords, test_sdrLegacyFilenameMapFullOnly

tests = functiontests(localfunctions);
end

function testMapIsWellFormed(testCase)
theMap = sdrLegacyFilenameMap();

collections = fieldnames(theMap.collections);
testCase.verifyEqual(sort(collections), ...
    sort({'cone_lattices'; 'cone_mosaics'; 'mrgc_lattices'; 'mrgc_on_circuits'}));

expectedCounts = struct('cone_lattices', 5, 'cone_mosaics', 27, ...
    'mrgc_lattices', 6, 'mrgc_on_circuits', 21);

for ii = 1:numel(collections)
    thisCollection = collections{ii};
    records = theMap.collections.(thisCollection).records;
    testCase.verifyEqual(numel(records), expectedCounts.(thisCollection), ...
        sprintf('Wrong record count for %s.', thisCollection));

    % Paths and checksums identify an asset, so both must be unique.
    testCase.verifyEqual(numel(unique({records.sdr_path})), numel(records), ...
        sprintf('Duplicate sdr_path in %s.', thisCollection));
    testCase.verifyEqual(numel(unique({records.sha256})), numel(records), ...
        sprintf('Duplicate sha256 in %s.', thisCollection));

    % A recorded legacy name must be unique too, or it cannot be a key.
    legacyNames = {records.legacy_filename};
    legacyNames = legacyNames(~cellfun(@isempty, legacyNames));
    testCase.verifyEqual(numel(unique(legacyNames)), numel(legacyNames), ...
        sprintf('Duplicate legacy_filename in %s.', thisCollection));

    for jj = 1:numel(records)
        testCase.verifyEqual(numel(records(jj).sha256), 64);
        testCase.verifyGreaterThan(records(jj).bytes, 0);
    end
end
end

function testEveryRemainingBundledFileMatchesItsRecord(testCase)
% Parity check. Run this before deleting any bundled asset.
theMap = sdrLegacyFilenameMap();
collections = fieldnames(theMap.collections);

checkedCount = 0;
for ii = 1:numel(collections)
    records = theMap.collections.(collections{ii}).records;
    for jj = 1:numel(records)
        thisRecord = records(jj);
        if isempty(thisRecord.legacy_repository_path), continue; end

        bundledFile = fullfile(isetbioRootPath, thisRecord.legacy_repository_path);
        if ~isfile(bundledFile), continue; end   % already migrated to the SDR

        details = dir(bundledFile);
        testCase.verifyEqual(details.bytes, thisRecord.bytes, ...
            sprintf('Byte count differs for %s.', thisRecord.legacy_repository_path));
        testCase.verifyEqual(lower(ieHash(bundledFile, 'file', 'SHA-256', 'hex')), ...
            lower(thisRecord.sha256), ...
            sprintf('SHA-256 differs for %s.', thisRecord.legacy_repository_path));
        checkedCount = checkedCount + 1;
    end
end

fprintf('Parity checked %d bundled files still present in the repository.\n', ...
    checkedCount);
end

function testUnresolvedCircuitsAreDeclared(testCase)
% Five deposited circuits never had a repository copy, so their legacy name
% could not be recovered by checksum. Each must still carry a candidate and
% an explicit confidence, so no caller mistakes a guess for a fact.
records = sdrLegacyFilenameMap('mrgc_on_circuits');

unresolved = records(cellfun(@isempty, {records.legacy_filename}));
testCase.verifyEqual(numel(unresolved), 5);

for ii = 1:numel(unresolved)
    testCase.verifyNotEmpty(unresolved(ii).candidate_legacy_filename);
    testCase.verifyTrue(ismember(unresolved(ii).candidate_confidence, ...
        {'high', 'unverified'}));
    testCase.verifyEqual(unresolved(ii).provenance, 'no-repository-copy');
end
end
