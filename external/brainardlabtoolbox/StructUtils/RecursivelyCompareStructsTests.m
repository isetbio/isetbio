function RecursivelyCompareStructsTests()
% RecursivelyCompareStructsTests
%
% Some basic unit tests for RecursivelyCompareStructs

% 2/7/17  NPC  Wrote it.

    %% Unit Test 1. Simplest case;

    s1 = struct('a', 1);
    s2 = struct('a', 1);
    result = RecursivelyCompareStructs('s1', s1, 's2', s2);

    s1 = struct('a', 1);
    s2 = struct('a', 1.02);
    result = RecursivelyCompareStructs('s1', s1, 's2', s2);

    %% More complicated examples
    % Make a filterBank struct with lots of nested stuff in it
    spatialFilter = struct(...
        'name', 'filterName', ...
        'center', 0.0, ...
        'sigma', 0.2, ...
        'phase', pi/3 ...
    );

    s1 = spatialFilter;
    s1.name = 'cosPhase';
    s1.phase = 0;

    s2 = spatialFilter;
    s2.name = 'sinPhase';
    s2.phase = pi/2;

    timeSamples = linspace(-1,1,100);
    t1 = struct(...
        'tau', 0.3 ...
    );

    % A struct array
    filterBank.spatialFilters = [s1, s2];
    
    % A cell array
    filterBank.temporalFilters = {'monophasic', timeSamples, t1, rand(1,10)};

    % A logical array
    filterBank.highRes = [true false];

    % A numeric array
    filterBank.spatialSamples = (-0.2:0.2:64);

    % A character array
    filterBank.name = 'Gabor';

    filterBank.otherInfo{1} = 'original';
    filterBank.otherInfo{2} = struct(...
            'isQuadrature', true, ...
            'coefficients', [1 0.2 nan eps 1e-13], ...
            'moreInfo', struct(...
                'useIsetbio', true, ...
                'subCoeffs', rand(1,3)));
    filterBank.otherInfo{3} = 1.0;


    %% Make  a perfect copy of filterBank
    filterBankCopy = filterBank;

    %% Set default toleance and empty customTolerances
    defaultTolerance = 1e-10;

    %% Unit Test 2 - identical structs
    RecursivelyCompareStructs('filterBank', filterBank, 'filterBankCopy', filterBankCopy, ...
           'defaultTolerance', defaultTolerance);
    
    %% Unit Test 3 - difference in a cell arrays
    filterBankCopy.otherInfo{1} = 'copy';
    filterBank.otherInfo{2}.moreInfo.useIsetbio = false;
    RecursivelyCompareStructs('filterBank', filterBank, 'filterBankCopy', filterBankCopy, ...
        'defaultTolerance', defaultTolerance);
    
    %% Unit Test 4- more differences in cell arrays
    filterBank.otherInfo{2}.moreInfo.subCoeffs(1) = filterBank.otherInfo{2}.moreInfo.subCoeffs(1) + 1e-5;
    RecursivelyCompareStructs('filterBank', filterBank, 'filterBankCopy', filterBankCopy, ...
        'defaultTolerance', defaultTolerance);
    
    %% Unit Test 5 - Apply custom tolerances
    customTolerances = containers.Map();
    customTolerances('otherInfo{2}.moreInfo.subCoeffs') =  1e-4;
    RecursivelyCompareStructs('filterBank', filterBank, 'filterBankCopy', filterBankCopy, ...
        'defaultTolerance', defaultTolerance, ...
        'customTolerances', customTolerances);
    
    %% Unit Test 6 - More custom tolerances
    customTolerances('temporalFilters{4}') = 2;
    filterBank.temporalFilters{4} = rand(1,10);
    RecursivelyCompareStructs('filterBank', filterBank, 'filterBankCopy', filterBankCopy, ...
        'defaultTolerance', defaultTolerance, ...
        'customTolerances', customTolerances);
end

