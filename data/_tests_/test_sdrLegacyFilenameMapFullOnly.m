function tests = test_sdrLegacyFilenameMapFullOnly()
% Live check that the legacy-filename map agrees with the published deposit
%
% Reads only the manifests, not the assets, so it is fast even though it
% touches the network. The FullOnly name keeps it out of the offline suite.
%
% See also
%   sdrLegacyFilenameMap, sdrMosaicRecords, test_sdrLegacyFilenameMap

tests = functiontests(localfunctions);
end

function testMapMatchesPublishedManifests(testCase)
theMap = sdrLegacyFilenameMap();
collections = fieldnames(theMap.collections);

for ii = 1:numel(collections)
    thisCollection = collections{ii};
    mapRecords = theMap.collections.(thisCollection).records;
    sdrRecords = sdrMosaicRecords(thisCollection);

    testCase.verifyEqual(numel(mapRecords), numel(sdrRecords), ...
        sprintf('%s record count differs from the deposit.', thisCollection));

    for jj = 1:numel(mapRecords)
        thisRecord = mapRecords(jj);
        sdrIndex = find(strcmp({sdrRecords.path}, thisRecord.sdr_path), 1);
        testCase.verifyNotEmpty(sdrIndex, ...
            sprintf('The deposit has no %s asset at %s.', ...
            thisCollection, thisRecord.sdr_path));

        testCase.verifyEqual(sdrRecords(sdrIndex).bytes, thisRecord.bytes);
        testCase.verifyEqual(lower(sdrRecords(sdrIndex).sha256), ...
            lower(thisRecord.sha256));
    end
end
end

function testCircuitEccentricitySignLossIsStillPresent(testCase)
% Guard on the known manifest defect. When a later deposit version fixes
% the sign, this test fails and the loader can stop using the legacy name
% as its lookup key.
sdrRecords = sdrMosaicRecords('mrgc_on_circuits');
negativeCount = sum(arrayfun(@(r) any(r.eccentricity_degrees < 0), sdrRecords));

testCase.verifyEqual(negativeCount, 0, ...
    ['The published mrgc_on_circuits manifest now records negative ', ...
    'eccentricities. Revisit sdrLegacyFilenameMap and the prebaked ', ...
    'circuit loader: selection by signed eccentricity may now be safe.']);
end
