function interpolatedValues = radiallyInterpolated2DMapFromMeridianValues(meridianValues, angles, interpolationMethod)
% Compute radially interpolated map from meridian map

   % Make sure all angles are > 0
   idx = find(angles<0);
   angles(idx) = angles(idx) + 360;
   
   % Add 360 deg angle and wrap around the values
   meridianAngles = RGCmodels.Watson.constants.indexedMeridianAngles;
   meridianAngles(5) = meridianAngles(1)+360;
   meridianValues(5,:) = meridianValues(1,:);
   
   % Do angular interpolation
   interpolatedValues = zeros(size(angles));
   parfor aa = 1:length(angles)
        interpolatedValues(aa) = interp1(meridianAngles, meridianValues(:,aa), angles(aa), interpolationMethod);
   end
   
end

