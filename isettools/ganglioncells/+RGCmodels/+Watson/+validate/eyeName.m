function eyeName(theEyeName)
% Validate the eye name

    m1 = RGCmodels.Watson.constants.leftEye;
    m2 = RGCmodels.Watson.constants.rightEye;
    assert(ismember(theEyeName, {m1, m2}), ...
         sprintf('eye name must be one of the following: ''%s'', ''%s''.', ...
         m1,m2));
    
end