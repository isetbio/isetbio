%% s_mosaicLibrary:  Standard sample mosaics
%
% Sometimes we just want a simple stored mosaic.  This code builds a
% library of a few mosaics with different sizes at different visual
% field positions
%
% We need to rerun this script every time we make a basic change to the
% @cMosaic class.  The saved mosaic may not know about the proper methods.
%
% See also
%   mosaicLoad

%% Sample sizes and positions

sizeDegs     = [0.5  0.5; 1.0 1.0; 2 2; 4 4;  8 8];  % multiple sizes
positionDegs = [0.0  0.0; 1.0 0.0; 2 0; 4 0; 10 0];  % Locations along the x-axis

dataDir = fullfile(isetRootPath,'data','cones');

%% Generate the mosaic

for ii=1:size(sizeDegs,1)
    for jj=1:size(positionDegs,1)

        cm = cMosaic(...
            'sizeDegs',     sizeDegs(ii,:), ...     % SIZE: 1.0 degs (x) 0.5 degs (y)
            'positionDegs', positionDegs(ii,:), ... % ECC: (0,0)
            'eccVaryingConeBlur', true ...
            );

        fname = mosaicName(sizeDegs,positionDegs);
        cm.save(fname,true);
        disp(fname)
    end
end

%% END