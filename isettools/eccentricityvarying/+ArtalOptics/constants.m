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
                     % Ranking strategy: 'resolution' (see analyzeArtalOptics.m)
                     ranking = [ ...
                          34    91    85    90    69    15 ...
                          33     7    46    60   130   128 ...
                          21    24    29    35    39    16 ...
                           2     8    73    37   102   105 ...
                          20    68     1    26    41    84 ...
                          80    45    95    62    13    82 ...
                          19   108    57    18    27     4 ...
                          52   119    28    48    97     5  ...
                          127    63    31    56    78    32];  
                         
                 case ArtalOptics.constants.rightEye
                     % Ranking strategy: 'resolution' (see analyzeArtalOptics.m)
                     ranking = [ ...
                         49    69   101    26    51    70  ...
                         56    65     7    15    81   102  ...
                         68     1    11   129   119    12  ...
                         88   110   111    82     8    22  ...
                         85    97   117    43     2    24  ...
                         45    41   125    99    37    80  ...
                         62    27   103   100    74     9  ...
                         18    34    31     6    84    14  ...
                         48    53    61    44    47    78  ...
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

