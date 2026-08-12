function tests = test_sdrMosaicCacheFullOnly()
% Live SDR regression test for the shared mosaic download and cache layer
%
% Exercises the cases every collection loader depends on: a first download,
% a cache hit, replacement of a corrupt cache entry, and a missing record.
% It also checks that the deposited lattices reproduce the golden crop
% values computed from the bundled lattices, which is the evidence needed
% before a bundled asset is deleted.
%
% Assets are removed from the local cache and downloaded again, so this test
% uses the network. The FullOnly name keeps it out of the offline suite.
%
% See also
%   sdrMosaicFetch, sdrMosaicRecords, sdrMosaicVerify, sdrLegacyFilenameMap

tests = functiontests(localfunctions);
end

function theFile = localFetchByLegacyName(collectionName, legacyName)
mapRecords = sdrLegacyFilenameMap(collectionName);
mapRecord = mapRecords(strcmp({mapRecords.legacy_filename}, legacyName));
assert(isscalar(mapRecord), 'No unique legacy record for %s.', legacyName);

sdrRecords = sdrMosaicRecords(collectionName);
sdrRecord = sdrRecords(strcmp({sdrRecords.path}, mapRecord.sdr_path));
theFile = sdrMosaicFetch(collectionName, sdrRecord);
end

function testColdFetchCacheHitAndCorruptReplacement(testCase)
% Use the smallest deposited cone mosaic so the test stays quick.
records = sdrMosaicRecords('cone_mosaics');
theRecord = records(strcmp({records.path}, ...
    'right/ecc-x0-y0__size-w0.5-h0.5.mat'));
testCase.verifyTrue(isscalar(theRecord));

cachedFile = fullfile(sdrMosaicCacheRoot(), 'cone_mosaics', 'right', ...
    'ecc-x0-y0__size-w0.5-h0.5.mat');

% Cold cache: the asset is downloaded and lands in the mirrored path.
if isfile(cachedFile), delete(cachedFile); end
fetchedFile = sdrMosaicFetch('cone_mosaics', theRecord);
testCase.verifyEqual(fetchedFile, cachedFile);
testCase.verifyTrue(sdrMosaicVerify(cachedFile, theRecord));

% Cache hit: same file, still verified, no temporary directories left over.
secondFile = sdrMosaicFetch('cone_mosaics', theRecord);
testCase.verifyEqual(secondFile, cachedFile);
leftovers = dir(fullfile(sdrMosaicCacheRoot(), 'tp*'));
testCase.verifyEmpty(leftovers, 'A temporary download directory was left behind.');

% Corrupt cache: a damaged file must be replaced, not returned.
fid = fopen(cachedFile, 'w');
fwrite(fid, zeros(1, 64, 'uint8'));
fclose(fid);
testCase.verifyFalse(sdrMosaicVerify(cachedFile, theRecord));

repairedFile = sdrMosaicFetch('cone_mosaics', theRecord);
testCase.verifyEqual(repairedFile, cachedFile);
testCase.verifyTrue(sdrMosaicVerify(cachedFile, theRecord));
end

function testMissingRecordsRaiseClearErrors(testCase)
testCase.verifyError(@() sdrMosaicRecords('no_such_collection'), ...
    'sdrMosaicRecords:UnknownCollection');
testCase.verifyError(@() mosaicLoad([99 99], [0 0]), ...
    'mosaicLoad:UnknownMosaic');
testCase.verifyError(@() retinalattice.import.sourceLatticeFile(...
    11, 'cones', 'left eye'), 'sourceLatticeFile:UnknownLattice');
end

function testDepositedLatticesReproduceGoldenCrops(testCase)
% The same golden values as test_mosaicLoadingGoldenValues, but computed
% from the deposited lattices rather than whatever is in the local gallery.
eccMicrons = 1000 * RGCmodels.Watson.convert.rhoDegsToMMs([0 0]);
sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(...
    [0.5 0.5], 0);

coneFile = localFetchByLegacyName('cone_lattices', ...
    'right_eye_cones_58deg_mosaic.mat');
coneData = load(coneFile, 'rfPositions');
conePositions = double(retinalattice.compute.croppedPositions(...
    -coneData.rfPositions, eccMicrons, sizeMicrons));

testCase.verifyEqual(size(conePositions), [3525 2]);
testCase.verifyEqual(conePositions([1 end], :), ...
    [66.8192901611328 55.5558128356934; ...
    -66.8488845825195 45.7090034484863], 'AbsTol', 1e-10);

mrgcFile = localFetchByLegacyName('mrgc_lattices', ...
    'right_eye_midget_ganglion_cells_64deg_mosaic.mat');
mrgcData = load(mrgcFile, 'rfPositions');
mrgcPositions = double(retinalattice.compute.croppedPositions(...
    -mrgcData.rfPositions, eccMicrons, sizeMicrons));

testCase.verifyEqual(size(mrgcPositions), [3265 2]);
testCase.verifyEqual(mrgcPositions([1 end], :), ...
    [66.8051986694336 39.4965705871582; ...
    -66.0864868164062 -64.4742279052734], 'AbsTol', 1e-10);
end

function testSmallestDepositedCircuitLoads(testCase)
% The 0.5-degree circuit is 2.7 MB, so the whole circuit download, checksum,
% and load path is exercised cheaply.
records = sdrMosaicRecords('mrgc_on_circuits');
theRecord = records(strcmp({records.path}, ...
    ['right/ecc-x2-y0__size-w0.5-h0.5__optics-polans2015-2__', ...
    'surround-packer-dacey-2002-h1-free-low.mat']));

theFile = sdrMosaicFetch('mrgc_on_circuits', theRecord);
loadedData = load(theFile, 'theMRGCMosaic');

testCase.verifyEqual(loadedData.theMRGCMosaic.rgcsNum, 557);
testCase.verifyEqual(loadedData.theMRGCMosaic.sizeDegs, [0.5 0.5], 'AbsTol', 1e-6);

% The manifest records this circuit at x = +2 degrees, and the asset agrees.
% The circuit at x = -2 degrees carries the same manifest eccentricity, which
% is why loaders select circuits by legacy filename. See sdrLegacyFilenameMap.
testCase.verifyEqual(loadedData.theMRGCMosaic.eccentricityDegs, [2 0], ...
    'AbsTol', 1e-6);
end
