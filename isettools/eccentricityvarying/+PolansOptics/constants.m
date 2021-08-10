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
       
        % Return subject ranking according to the PSF resolution at the
        % fovea. In this ranking, the higest resolving subjects appear
        % first. Ranking is computed by a call to analyzePolansOptics(true)
        % found in isetbio/calculators/opticsAssessment
        function  ranking = subjectRanking
            % Ranking strategy: 'resolution' (see analyzePolansOptics.m)
            ranking = [...
                6      9     2  ...
                10     8     4  ...
                3      5     1  ...
                7];
        end
        
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
        
        function flag = subjectRequiresCentralRefractionCorrection(subjectID) 
            subjectsRequiringCentralRefractionCorrection = setdiff(1:10, [1 2 4 7]);
            flag = ismember(subjectID, subjectsRequiringCentralRefractionCorrection);
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
            x = [];
            y = [];
            v = [];
            xq = [];
            yq = [];
            for vEccIndex = 1:size(zMap,1)
                for hEccIndex = 1:size(zMap,2)
                    % Outliers: exclude raw data and interpolate from neighboring data
                    
                   % ------ Subject 1 outliers ------ 
                    inBadPointsSetOfSubject1 = ...
                        (subjectIndex == 1) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 0) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-15 -14 -13 -12 -11 -9 -8 -7 -6]));
                  
                    
                    % ------ Subject 2 outliers ------ 
                    inBadPointsSet1OfSubject2 = ...
                        (subjectIndex == 2) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 0) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-7 -6 -5 -4  -2 -1 0 1]));
                    inBadPointsSet2OfSubject2 = ...
                        (subjectIndex == 2) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -5) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-1 1 2 3]));
                    inBadPointsSetOfSubject2 = inBadPointsSet1OfSubject2 | inBadPointsSet2OfSubject2;
                    
                    
                    % ------ Subject 3 outliers ------ 
                    inBadPointsSetOfSubject3= ...
                        (subjectIndex == 3) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 0) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-1]));
                    
                    
                    % ------ Subject 5 outliers ------ 
                    inBadPointsSet1OfSubject5= ...
                        (subjectIndex == 5) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -5) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-7 -6 -5]));
                    
                    inBadPointsSet2OfSubject5= ...
                        (subjectIndex == 5) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -10) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-15]));
                    
                    inBadPointsSet3OfSubject5= ...
                        (subjectIndex == 5) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 0) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-9 -10]));
                    
                    inBadPointsSetOfSubject5 = inBadPointsSet1OfSubject5 | inBadPointsSet2OfSubject5 | inBadPointsSet3OfSubject5;
                   
                    
                    % ------ Subject 6 outliers ------ 
                    inBadPointsSet1OfSubject6 = ...
                        (subjectIndex == 6) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 0) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-9 -8 -7 -6]));
                    
                    inBadPointsSet2OfSubject6 = ...
                        (subjectIndex == 6) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -5) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-9 -7 -6 7]));
                    
                    inBadPointsSet3OfSubject6 = ...
                        (subjectIndex == 6) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -10) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-9 -7 -6 7]));
                 
                     inBadPointsSet4OfSubject6 = ...
                        (subjectIndex == 6) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 5) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-10 -12]));
                    inBadPointsSetOfSubject6 = inBadPointsSet1OfSubject6 | inBadPointsSet2OfSubject6 | inBadPointsSet3OfSubject6 | inBadPointsSet4OfSubject6;
                    
                    
                    % ------ Subject 9 outliers ------ 
                    inBadPointsSet1OfSubject9 = ...
                        (subjectIndex == 9) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -15) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [-5]));
                    
                    inBadPointsSet2OfSubject9 = ...
                        (subjectIndex == 9) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == 5) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [10]));
                    inBadPointsSetOfSubject9 = inBadPointsSet1OfSubject9 | inBadPointsSet2OfSubject9;
                    
                    
                    % ------ Subject 10 outliers ------ 
                    inBadPointsSetOfSubject10 = ...
                        (subjectIndex == 10) && ...
                        (PolansOptics.constants.measurementVerticalEccentricities(vEccIndex) == -10) && ...
                        (ismember(PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex), [0 1 2 3 4]));
                    
                    
                    xyIndex = xyIndex+1;
                    theCoeffs = zCoeffs(xyIndex,:);
                    zMap(vEccIndex, hEccIndex,:) = theCoeffs;
                    % Keep the good (x,y) points
                    if ((~all(theCoeffs == 0)) && (~inBadPointsSetOfSubject1) && (~inBadPointsSetOfSubject2) && (~inBadPointsSetOfSubject3) && ...
                            (~inBadPointsSetOfSubject5) && (~inBadPointsSetOfSubject6) && (~inBadPointsSetOfSubject9) && ...
                            (~inBadPointsSetOfSubject10))
                        x = cat(2, x, PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex));
                        y = cat(2, y, PolansOptics.constants.measurementVerticalEccentricities(vEccIndex));
                        v = cat(2, v, theCoeffs');
                    end
                    xq = cat(2, xq, PolansOptics.constants.measurementHorizontalEccentricities(hEccIndex));
                    yq = cat(2, yq, PolansOptics.constants.measurementVerticalEccentricities(vEccIndex));
                end
            end
            
            
            for z = 1:size(zMap,3)
                % Fill in all-zero entries
                F = scatteredInterpolant(x',y',(v(z,:))');
                tmp = F(xq,yq);
                zMapFilled(:,:,z) = (reshape(tmp, [size(zMap,2), size(zMap,1)]))';
            end

            % Override with filled zmap
            zMap = zMapFilled;
        end
        

    end % Static methods
    
end
