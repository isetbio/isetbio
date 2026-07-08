function meridianName(theMeridianName)
% Validate the meridian name

    m1 = RGCmodels.Watson.constants.nasalMeridian;
    m2 = RGCmodels.Watson.constants.temporalMeridian;
    m3 = RGCmodels.Watson.constants.superiorMeridian;
    m4 = RGCmodels.Watson.constants.inferiorMeridian;
    assert(ismember(theMeridianName, {m1, m2, m3, m4}), ...
         sprintf('meridian name must be one of the following: ''%s'', ''%s'', ''%s'', ''%s''.', ...
         m1,m2,m3,m4));
    
end

