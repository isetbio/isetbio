function [coneType, filters, filterNames] = mouseConeMosaic(sz,fHeight, rgbDensities, sensor)
%Create a spatial representation of the mouse cone mosaic
% This divides the sensor into strips : M, UV+M and UV.
% WARNING : the division is defined as a fraction of the sensor size, so
% changing the fov will change the retina's mosaic!! fov=30 is good.
%
% Drafted by EC

if notDefined('sz'), error('Array size must be defined'); end

if notDefined('sensor')
    % load filters data
    mouseFilters = load('~/psych221/mouseColorFilters.mat');
else
    mouseFilters.data = sensorGet(sensor, 'filterSpectra');
    mouseFilters.filterNames = sensorGet(sensor, 'filterNames');
    mouseFilters.wavelength = sensorGet(sensor, 'wave');
end

if notDefined('fHeight')
    fil =  load('~/psych221/mouseColorFilters.mat');
    fHeight = fil.fHeight; % default : [-0.5, -0.1, 0.1, 0.5]
    clear fil;
end
if rgbDensities(1)==0 && rgbDensities(2)==0
    error('no mouse cones at all!')
elseif rgbDensities(1)==0
    % monochrome with only UV cones
    filters = mouseFilters.data;
    filterNames = mouseFilters.filterNames;
    coneType = 2*ones(sz);
elseif rgbDensities(2)==0
    % monochrome with only M cones
    filters = mouseFilters.data;
    filterNames = mouseFilters.filterNames;
    coneType = ones(sz);
else
    % normal polychrome
    % Place M cones on the top, mixed cones in the middle and UV cones at the
    % bottom. TODO : orient the mosaic with an angle.
    fHeight = fHeight/max(fHeight)/2 + 0.5; % convert to [0;1]
    fHeight = floor(fHeight*sz(1)); % convert to height in image coordinates
    
    % number of lines of mixed cones :
    numMixed = fHeight(3)-fHeight(2);
    
    % set the map of cone types
    coneType = zeros(sz);
    coneType( 1 : fHeight(2), :) = 1; % M-cones
    for i=1:numMixed
        coneType( (fHeight(2)+i) , :) = 1+i; % mixed cones; TODO interpolation
    end
    coneType( (fHeight(3)+1) : sz(1)     , :) = 2 + numMixed; % UV-cones
    
    % Create the filters
    filters = zeros(length(mouseFilters.wavelength), numMixed+2);
    filters(:,1) = mouseFilters.data(:,1);
    filters(:,end) = mouseFilters.data(:,2);
    interpWeights = (0:numMixed+1)/(numMixed+1);
    interpWeights = interpWeights(2:end-1);
    filters(:,2:end-1) = filters(:,end) * interpWeights + filters(:,1) * (1-interpWeights);
    
    % interpolation : the band in the middle of the retina is a mix of UV+M
    % pigments. We use linear combinations of the UV and M spectra. This
    % creates (a possibly large number of) new cone types.
    % The cone types / color filters are numbered in order, from top to bottom,
    % which coincides with the "filters" array.
    
    % Create the filter names
    numFilters = numMixed+2;
    filterNames = cell(1,numFilters);
    filterNames{1} = mouseFilters.filterNames{1};
    for i=2:(numFilters-1)
        filterNames{i} = sprintf('mouse_UV+M_cone_%f', interpWeights(i-1));
    end
    filterNames{numFilters} = mouseFilters.filterNames{2};
    
end

end