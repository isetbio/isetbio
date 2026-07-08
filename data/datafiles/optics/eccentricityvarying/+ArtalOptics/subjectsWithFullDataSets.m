function [targetEyeSubjectPool, binocularSubjectPool] = subjectsWithFullDataSets(whichEye)

   targetEcc = 0;
   pupilDiamMM = 3;
   wavelengthsListToCompute = 550;
   micronsPerDegree = 300;
   
   targetEyeSubjectPool= [];
   binocularSubjectPool = [];

   for subjectID = 1:130
       theOI = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                    whichEye, targetEcc, pupilDiamMM, ...
                    wavelengthsListToCompute, micronsPerDegree);
       if (contains(whichEye, 'right')) ...
          theOtherOI = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                    'left eye', targetEcc, pupilDiamMM, ...
                    wavelengthsListToCompute, micronsPerDegree);
       else
          theOtherOI = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                    'right eye', targetEcc, pupilDiamMM, ...
                    wavelengthsListToCompute, micronsPerDegree);
       end

       if (isempty(theOI))
           %fprintf(2,'Nan or zeros for subject %d and eye %s\n', ...
           %    subjectID, whichEye);
       else
           targetEyeSubjectPool = cat(2, targetEyeSubjectPool, subjectID);
           if (~isempty(theOtherOI))
               binocularSubjectPool = cat(2, binocularSubjectPool, subjectID);
           end
       end
       
   end
   
end