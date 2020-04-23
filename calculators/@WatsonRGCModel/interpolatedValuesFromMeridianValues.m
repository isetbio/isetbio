function val = interpolatedValuesFromMeridianValues(obj, meridianValues, requestedAngles)
% Angular interpolation from meridian values to full 360

   % Make sure all angles are > 0
   idx = find(requestedAngles<0);
   requestedAngles(idx) = requestedAngles(idx) + 360;
   
   % Add 360 point and wrap around the values
   meridianAngles = obj.enumeratedMeridianAngles;
   meridianAngles(5) = meridianAngles(1)+360;
   meridianValues(5,:) = meridianValues(1,:);
   
   % Do the angular interpolation
   val = zeros(size(requestedAngles));
   method = 'linear'; %'makima'; % 'linear'; % 'spline'; % 'linear'
   parfor aa = 1:length(requestedAngles)
        val(aa) = interp1(meridianAngles, meridianValues(:,aa), requestedAngles(aa), method);
   end
   
end
