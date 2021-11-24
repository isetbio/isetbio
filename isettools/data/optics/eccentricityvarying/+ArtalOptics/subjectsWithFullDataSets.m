function goodSubjectPool = subjectsWithFullDataSets(whichEye)

   targetEcc = 0;
   pupilDiamMM = 3;
   wavelengthsListToCompute = 550;
   micronsPerDegree = 300;
   
   goodSubjectPool = [];
   for subjectID = 1:130
       theOI = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                    whichEye, targetEcc, pupilDiamMM, ...
                    wavelengthsListToCompute, micronsPerDegree);
       if (isempty(theOI))
           %fprintf(2,'Nan or zeros for subject %d and eye %s\n', ...
           %    subjectID, whichEye);
       else
           goodSubjectPool = cat(2, goodSubjectPool, subjectID);
       end
       
   end
   
end