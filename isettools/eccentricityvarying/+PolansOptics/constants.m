classdef constants
% Constants from Polans study

    properties (Constant)
        % Pupil diameter for all measurements
        measurementPupilDiamMM = 4.0;
          
        % Eccentricities at which measurements were made
        measurementHorizontalEccentricities = 40:-1:-40;
        measurementVerticalEccentricities = 25:-5:-25;
    end
    
    methods (Static)
        function [zMap, zCoeffIndices] = ZernikeCoefficientsMap(subjectIndex)
            % Import raw data
            allData = rawDataReadData('zCoefsPolans2015', 'datatype', 'isetbiomatfileonpath');
            allData = allData.data;
    
            % Reported Z-coeffs are Z3-Z20
            zCoeffIndices = 3:size(allData,3);
            zCoeffs = squeeze(allData(subjectIndex , :, zCoeffIndices));
            
            % Assemble xy map of Zernike coefficients
            zMap = zeros(...
                numel(PolansOptics.constants.measurementVerticalEccentricities), ...
                numel(PolansOptics.constants.measurementHorizontalEccentricities), ...
                numel(zCoeffIndices));
            
            xyIndex = 0;
            for vEccIndex = 1:size(zMap,1)
                for hEccIndex = 1:size(zMap,2)
                    xyIndex = xyIndex+1;
                    zMap(vEccIndex, hEccIndex,:) = zCoeffs(xyIndex,:);
                end
            end
        end
    end % Static methods
    
end
