classdef constants
% Constants

    properties (Constant)

        LcenterColor = [240 185 192]/255;
        McenterColor = [100 210 220]/255;
        ScenterColor = [200 100 250]/255;
        achromaticColor = [0.8 0.8 0.8];

        validPackerDacey2002H1CellIndices = [1 2 3 4];

        PackerDaceyParamsForH1CellIndices = struct(...
                    'RnarrowToRwideRatios',  [152/515 170/718 115/902 221/1035], ...
                    'VnarrowToVwideRatios',  [1.0     0.8     0.3     0.2], ...
                    'source', '4 cells in Figure 6 of Packer & Dacey (2002)' ...
                    );

    end
end
