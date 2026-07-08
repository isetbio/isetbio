classdef constants
% Constants from the Thibos study

     properties (Constant)
        % Label for left eye
        leftEye = 'left eye';
        
        % Label for right eye
        rightEye = 'right eye';

        % Default measurement pupil diameter
        measurementPupilDiamMM = 7.5;
        
        % Available measurement pupil diameters
        availableMeasurementPupilDiamsMM = [3.0 4.5 6.0 7.5];

        % Data exist for both left and right eyes
        measurementEyeOrder = {'right eye', 'left eye'};
     
     end


     methods (Static)
         
        % Return subject ranking according to the PSF resolution at the
        % fovea. In this ranking, the higest resolving subjects appear
        % first. Ranking is computed by a call to analyzeThibosOptics(true)
        % found in isetbio/calculators/opticsAssessment
        function ranking = subjectRanking(whichEye)
             switch (whichEye)
                 case ThibosOptics.constants.leftEye
                    % Ranking strategy: 'resolution' (see analyzeThibosOptics.m)
                     ranking = [ ...
                          27    66    67     2    26    55     1  ...
                          62    49    43    34    25    28    24  ...
                          17    39    21    52    10    16    63  ...
                          33    38    59    22    42    53     5  ...
                          11    14     4    35    19    56     6  ...
                          37    40    60    29    30    18    44  ...
                          32    20    31    36    48    65     8  ...
                          50    15    46    45    58    41     7  ...
                          13    61    47     9    54    23    64  ...
                           3    57    12    51    69    68    70];


                case ThibosOptics.constants.rightEye
                     % Ranking strategy: 'resolution' (see analyzeThibosOptics.m)
                     ranking = [ ...
                          26    62    43    66     5    31    47  ...
                          14    38    34    48    21    40    50  ...
                           3    58    35    17    16    59    27  ...
                          65    42    41    63    49     8    36  ...
                           1    24    33    55    29    22     2  ...
                          10    28    11    53    18    19    57  ...
                          52    32    25     4     7    54     6  ...
                          60    45    67    37    12    39    56  ...
                          30    51    15    46    20    44     9  ...
                          23    13    61    64    68    69    70];
                         
             end    
        end

        function flag = subjectRequiresCentralRefractionCorrection(whichEye, subjectID)
            % Not implemented yet
            flag = false;
        end

        function [zMap, zCoeffIndices] = ZernikeCoefficientsMap(subjectIndex, whichEye)
        end

     end

end
