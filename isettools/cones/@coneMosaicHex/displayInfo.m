function displayInfo(obj, varargin)
% Print various infos about the cone mosaic
%
% Syntax:
%    displayInfo(obj)
%
% Description:
%    Print various infos about the cone mosaic
%
% Inputs:
%    obj - The cone mosaic hex object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    plotApertureStats         - Boolean, whether the plot a histogram of
%                                aperture diameter and light collecting area

% History:
%    xx/xx/15  NPC  ISETBIO TEAM, 2015
%    02/16/18  jnm  Formatting
%    10/12/18  NPC  Now reporting & plotting stats for inner segment aperture and area

% parse input
p = inputParser;
p.addParameter('plotApertureStats', false, @islogical);
p.parse(varargin{:});

if (isempty(obj.apertureStats))
    plotApertureStats = false;
    obj.computeApertureStats(plotApertureStats);
end


[~,~,~,maximumConeDensity, minimumConeSeparationMicrons] = ...
    obj.computeDensityMap('from mosaic');
        
fprintf('\nMosaic info:\n');
fprintf('%53s %2.1f (w) x %2.1f (h)\n', 'Size (microns):', ...
    obj.width * 1e6, obj.height * 1e6);
fprintf('%53s %2.2f (w) x %2.2f (h)\n', 'FOV (deg):', ...
    obj.fov(1), obj.fov(2));
fprintf('%53s %2.3f (aperture) x %2.3f (geometric)\n', ...
    'Retinal coverages:', obj.innerSegmentCoverage, obj.coverage);
fprintf('%53s %0.0f\n', 'Grid resampling factor:', obj.resamplingFactor);

if (obj.eccBasedConeDensity == false)
    fprintf('%53s %2.2f (w) x %2.2f (h), diameter: %2.2f\n', ...
        'Cone geometric aperture (microns):', ...
        obj.pigment.width * 1e6, obj.pigment.height * 1e6, ...
        diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.width * 1e6));
    fprintf('%53s %2.2f (w) x %2.2f (h), diameter: %2.2f\n', ...
        'Cone light colleting aperture (microns):', ...
        obj.pigment.pdWidth * 1e6, obj.pigment.pdHeight * 1e6, ...
        diameterForCircularApertureFromWidthForSquareAperture(obj.pigment.pdWidth * 1e6));
else
    fprintf('%53s Min=%2.4f, Mean=%2.4f, Median=%2.4f, Max=%2.4f \n', 'Inner segment diameter (microns): ', ...
        obj.apertureStats.rangeDiameterMicrons(1), obj.apertureStats.meanDiameterMicrons, obj.apertureStats.medianDiameterMicrons, obj.apertureStats.rangeDiameterMicrons(2));
    fprintf('%53s Min=%2.4f, Mean=%2.4f, Median=%2.4f, Max=%2.4f \n', 'Inner segment area (microns^2): ', ...
        obj.apertureStats.rangeLightCollectingArea(1), obj.apertureStats.meanLightCollectingArea, obj.apertureStats.medianLightCollectingArea, obj.apertureStats.rangeLightCollectingArea(2));
end


fprintf('%53s %2.4f \n', 'Cone geometric area (microns^2):', ...
    obj.pigment.area * 1e12);
fprintf('%53s %2.4f\n', 'Cone light colleting area (microns^2):', ...
    obj.pigment.pdArea * 1e12);

%fprintf('%53s %2.3f \n', 'Cone coverage :', obj.coverage);
%fprintf('%53s %2.3f \n', 'Cone coverage (inner segments):', obj.innerSegmentCoverage);
fprintf('%53s %2.0f cols x %2.0f rows\n', 'Rectangular grid:', ...
    size(obj.patternOriginatingRectGrid, 2), ...
    size(obj.patternOriginatingRectGrid, 1));
fprintf('%53s %2.0f cols x %2.0f rows\n', 'Resampled grid:', ...
    obj.cols, obj.rows);
fprintf('%53s %d\n', 'Total cones:', numel(obj.pattern));
totalConesNum = numel(find(obj.pattern > 1));
LconesNum = numel(find(obj.pattern == 2));
MconesNum = numel(find(obj.pattern == 3));
SconesNum = numel(find(obj.pattern == 4));
fprintf('%53s %d\n', 'Active cones:' , totalConesNum);
fprintf('%53s %d (%2.3f%%)\n', 'L- cones:' , LconesNum, 100*LconesNum/totalConesNum);
fprintf('%53s %d (%2.3f%%)\n', 'M- cones:' , MconesNum, 100*MconesNum/totalConesNum);
fprintf('%53s %d (%2.3f%%)\n', 'S- cones:' , SconesNum, 100*SconesNum/totalConesNum);
fprintf('%53s %2.1f cones/mm^2\n', 'Cone density (all cones):', ...
    numel(obj.pattern) / (obj.width * obj.height * 1e6));
fprintf('%53s %2.1f cones/mm^2\n', 'Mean cone density (active cones):', ...
    numel(find(obj.pattern > 1)) / (obj.width * obj.height * 1e6));

fprintf('%53s %2.1f cones/mm^2 at %2.2f,%2.2f um\n', 'Max cone density (based on 5 neighboring cones):', ...
maximumConeDensity.value, maximumConeDensity.position(1), maximumConeDensity.position(2));
fprintf('%53s %2.2f microns at %2.2f,%2.2f um\n', 'Min cone separation (based on 5 neighboring cones):', ...
minimumConeSeparationMicrons.value, minimumConeSeparationMicrons.position(1), minimumConeSeparationMicrons.position(2));
fprintf('%53s %d\n', 'Ecc-based cone efficiency:', obj.eccBasedConeQuantalEfficiency);
fprintf('%53s %d\n\n', 'Ecc-based macular pigment density:', obj.eccBasedMacularPigment);
end
