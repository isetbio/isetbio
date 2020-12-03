classdef constants
% Constants from Polans study

    properties (Constant)
        % Label for left eye
        leftEye = 'left eye';
        
        % Label for right eye
        rightEye = 'right eye';
        
        % Label for nasal meridian
        nasalMeridian    = 'nasal meridian';
        
        % Label for temporal meridian
        temporalMeridian = 'temporal meridian';
        
        % Label for inferior meridian
        inferiorMeridian = 'inferior meridian';
        
        % Label for superior meridian
        superiorMeridian = 'superior meridian';
        
        % Pupil diameter for all measurements
        measurementPupilDiamMM = 4.0;
        
        % All data are for the right eye
        measurementEye = 'right';
        
        % Eccentricities at which measurements were made       
        measurementHorizontalEccentricities = 40:-1:-40;
        measurementVerticalEccentricities = 25:-5:-25;
    end
   
    methods (Static)
        function retinalQuadrantName = retinalQuadrantFromEcc(ecc, whichEye)
            assert(numel(ecc) == 2, 'Eccentricity must be a 2-element vector (x,y)');
            assert(ismember(whichEye, {PolansOptics.constants.leftEye, PolansOptics.constants.leftEye}));
            switch (whichEye)
                case 'right'
                    if (ecc(1) < 0)
                        horizontalMeridianName = PolansOptics.constants.nasalMeridian;
                    elseif (ecc(1) == 0)
                        horizontalMeridianName = 'foveal pit';
                    else
                        horizontalMeridianName = PolansOptics.constants.temporalMeridian;
                    end
                case 'left'
                    if (ecc(1) < 0)
                        horizontalMeridianName = PolansOptics.constants.temporalMeridian;
                    elseif (ecc(1) == 0)
                        horizontalMeridianName = 'foveal pit';
                    else
                        horizontalMeridianName = PolansOptics.constants.nasalMeridian;
                    end
            end
            if (ecc(2) < 0)
                verticalMeridianName = PolansOptics.constants.inferiorMeridian;
            elseif (ecc(1) == 0)
                verticalMeridianName = 'foveal pit';
            else
                verticalMeridianName = PolansOptics.constants.superiorMeridian;
            end
            retinalQuadrantName{1} = horizontalMeridianName;
            retinalQuadrantName{2} = verticalMeridianName;
        end
        
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
