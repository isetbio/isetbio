function tests = test_dataPaths()
tests = functiontests(localfunctions);
end

function testCanonicalDataPath(testCase)
testFile = mfilename('fullpath');
expectedRoot = fileparts(fileparts(fileparts(testFile)));
expectedPath = fullfile(expectedRoot, 'data', 'datafiles');

testCase.verifyEqual(isetbioRootPath, expectedRoot);
testCase.verifyEqual(isetbioDataPath, expectedPath);
testCase.verifyTrue(isfolder(isetbioDataPath));
testCase.verifyFalse(isfolder(fullfile(isetbioRootPath, 'dataiset')));
end

function testRepresentativeBundledData(testCase)
relativeFiles = { ...
    fullfile('human', 'cones', 'coneDensityCurcio1990.mat'); ...
    fullfile('human', 'wvf', 'zCoefsJaekenArtal2012.mat'); ...
    fullfile('rgc', 'parasolData.mat'); ...
    fullfile('treeshrew', 'treeshrewLensAbsorbance.mat')};

for ii = 1:numel(relativeFiles)
    testCase.verifyTrue(isfile(fullfile(isetbioDataPath, relativeFiles{ii})), ...
        sprintf('Missing bundled data file: %s', relativeFiles{ii}));
end
end

function testDataReadersUseCanonicalTree(testCase)
mosaicFile = mosaicName([1 1], [0 0]);
testCase.verifyTrue(startsWith(mosaicFile, [isetbioDataPath filesep]));

% Source lattices now come from the isetbio-mosaics SDR deposit, so the
% gallery directories are empty in a fresh checkout. What must stay stable is
% where a locally generated lattice is written and looked for, and that the
% deposit still publishes the lattices the standard configurations name.
% Checking availability here rather than downloading keeps this test offline.
coneLatticeParams = retinalattice.configure(64, 'cones', 'right eye');
testCase.verifyEqual(coneLatticeParams.latticeGalleryDir, ...
    fullfile(isetbioDataPath, 'cones', 'lattices'));
testCase.verifyTrue(isfolder(coneLatticeParams.latticeGalleryDir));
verifyDepositPublishes(testCase, 'cone_lattices', ...
    coneLatticeParams.patchFinalPositionsSaveFileName);

mRGCLatticeParams = retinalattice.configure(60, ...
    'midget ganglion cells', 'right eye');
testCase.verifyEqual(mRGCLatticeParams.latticeGalleryDir, ...
    fullfile(isetbioDataPath, 'rgc', 'lattices'));
testCase.verifyTrue(isfolder(mRGCLatticeParams.latticeGalleryDir));
verifyDepositPublishes(testCase, 'mrgc_lattices', ...
    mRGCLatticeParams.patchFinalPositionsSaveFileName);

% The default mRGC synthesis configuration uses a 64-degree source lattice.
% Keep this check separate from the 60-degree compatibility asset above so a
% data-layout change cannot silently break standard circuit synthesis.
mRGCLattice64Params = retinalattice.configure(64, ...
    'midget ganglion cells', 'right eye');
testCase.verifyEqual(mRGCLattice64Params.latticeGalleryDir, ...
    fullfile(isetbioDataPath, 'rgc', 'lattices'));
verifyDepositPublishes(testCase, 'mrgc_lattices', ...
    mRGCLattice64Params.patchFinalPositionsSaveFileName);

coneDensity = rawDataReadData('coneDensityCurcio1990', ...
    'datatype', 'isetbiomatfileonpath');
testCase.verifyNotEmpty(fieldnames(coneDensity));

packageFunction = which('ArtalOptics.constants');
testCase.verifyTrue(startsWith(packageFunction, [isetbioDataPath filesep]));
end

function verifyDepositPublishes(testCase, collectionName, legacyFileName)
% The deposit must still publish the lattice a standard configuration names.
records = sdrLegacyFilenameMap(collectionName);
testCase.verifyTrue(any(strcmp({records.legacy_filename}, legacyFileName)), ...
    sprintf('The %s collection no longer publishes %s.', ...
    collectionName, legacyFileName));
end

function testBundledDataGoldenValues(testCase)
%% Spot-check representative data from both former data directories.

numericTolerance = 1e-12;

coneDensity = load(fullfile(isetbioDataPath, 'human', 'cones', ...
    'coneDensityCurcio1990.mat'));
testCase.verifySize(coneDensity.superior.density, [33 1]);
testCase.verifyEqual(coneDensity.superior.density([1 end]), ...
    [250000; 3239.4657119250924], 'AbsTol', numericTolerance);
testCase.verifyEqual(coneDensity.temporal.eccMM([1 end]), [0; 18], ...
    'AbsTol', numericTolerance);

wavefront = load(fullfile(isetbioDataPath, 'human', 'wvf', ...
    'zCoefsJaekenArtal2012.mat'));
testCase.verifySize(wavefront.data, [3901 84]);
testCase.verifyEqual(wavefront.data(1000, 20), -0.25988552258270098, ...
    'AbsTol', numericTolerance);
testCase.verifyEqual(wavefront.data(end, end), -0.0098923178097894901, ...
    'AbsTol', numericTolerance);

rgc = load(fullfile(isetbioDataPath, 'rgc', 'parasolData.mat'));
testCase.verifySize(rgc.parasolData, [187 2]);
testCase.verifyEqual(rgc.parasolData([1 end], :), ...
    [17.281865 293.435449; 3.031088 78.774617], ...
    'AbsTol', numericTolerance);

treeShrew = load(fullfile(isetbioDataPath, 'treeshrew', ...
    'treeshrewLensAbsorbance.mat'));
testCase.verifyEqual(treeShrew.wavelength([1 186 end]), [325 510 695]);
testCase.verifyEqual(treeShrew.data([1 186 end]), ...
    [0.7911817389734721 0.041072628226068955 0.038942918367977378], ...
    'AbsTol', numericTolerance);

end

function testISETCamDataRemainsSeparate(testCase)
relativeFiles = { ...
    fullfile('human', 'stockman.mat'); ...
    fullfile('human', 'luminosity.mat'); ...
    fullfile('human', 'XYZ.mat')};

for ii = 1:numel(relativeFiles)
    testCase.verifyTrue(isfile(fullfile(isetRootPath, 'data', relativeFiles{ii})), ...
        sprintf('Missing required ISETCam data file: %s', relativeFiles{ii}));
end
end
