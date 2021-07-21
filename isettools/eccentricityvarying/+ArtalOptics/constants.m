classdef constants
% Constants from the Jaeken and Artal, 2012. study

     properties (Constant)
        % Label for left eye
        leftEye = 'left eye';
        
        % Label for right eye
        rightEye = 'right eye';
        
        % Pupil diameter for all measurements
        measurementPupilDiamMM = 4.0;
        
        % Data exist for both left and right eyes
        measurementEyeOrder = {'right eye', 'left eye'};
        
        % Eccentricities at which measurements were made       
        measurementHorizontalEccentricities = -40:1:40;
        measurementVerticalEccentricities = 0;
     end
     
     
     methods (Static)
         
        % Return subject ranking according to the PSF resolution at the
        % fovea. In this ranking, the higest resolving subjects appear
        % first. Ranking is computed by a call to analyzeArtalOptics(true)
        % found in isetbio/calculators/opticsAssessment
         function ranking = subjectRanking(whichEye)
             switch (whichEye)
                 case ArtalOptics.constants.leftEye
                     ranking = [ ...
                         69     8    80    21    24    60   ...
                         73    85   102    37     1   128   ...
                         29    62   105    16    34   130   ...
                         90    46    45    68    33    95   ...
                         7     78    91     2    35    15   ...
                         26    28    39    20    63    84   ...
                         119    4    13    41    82    19   ...
                         48    97    18   108    57   127   ...
                         27    52     5    31    56    32];  
                         
                 case ArtalOptics.constants.rightEye
                     ranking = [ ...
                         49   102    62    68    81    69 ...
                         51   101    88   110    74     1 ...
                         26    11    70     8    56   129  ...
                         65   103     7   119    85    45  ...
                         15   100    82    80    34     2  ...
                         18    14    24    97   111   117  ...
                         22    37   125    12     9    61  ...
                         43    99    41    44    48    27  ...
                         6    53    84    31    47    78   ...
                         32];
             end
                     
         end
         
         function flag = subjectRequiresCentralRefractionCorrection(whichEye, subjectID) 
             
            lowResLeftEyeSubjects = [2     5    15    27    31    32    33    52    56    57    68    73    90   108   127   130];
            lowResRightEyeSubjects = [6     9    11    12    15    22    27    31    32    41    43    47    49    51    53    56    68    78    84    97    99   101   111   125];
            
            switch (whichEye)
                case 'left eye'
                    flag = ismember(subjectID, lowResLeftEyeSubjects);
                case 'right eye'
                    flag = ismember(subjectID, lowResRightEyeSubjects);
            end
            
            if (flag)
                fprintf('Subject %d from %s requires central refraction correction\n', subjectID, whichEye);
            end
            
         end
        
         function [zMap, zCoeffIndices] = ZernikeCoefficientsMap(subjectIndex, whichEye)
            % Import raw data
            allData = rawDataReadData('zCoefsJaekenArtal2012', 'datatype', 'isetbiomatfileonpath');
            allData = allData.data;
            
            % Horizontal eccentricities
            horzEccen = allData(1, 4:end);
            
            % Subjects num
            subjectsNum = sum(~isnan(unique(allData(:, 1))));
            
            % Reported Z-coeffs are Z0-Z14
            zCoeffIndices = 0:14;
            
            % Remove headers 
            allData = allData(2:end, 4:end);
            
            % Reshape the data: zernike x subject x eye (right, left) x eccentricity
            allData = reshape(allData, length(zCoeffIndices), subjectsNum, 2, numel(horzEccen));
            
            eyeIndex = strcmp(whichEye, ArtalOptics.constants.measurementEyeOrder);
            zMap = (squeeze(allData(:, subjectIndex, eyeIndex,:)))';
          
         end
         
     end
     
end

