function [primary] = SettingsToPrimary(calOrCalStruct,settings)
% [primary] = SettingsToPrimary(calOrCalStruct,settings)
% 
% Convert from device settings coordinates to
% primary coordinates by inverting
% the gamma correction.
%
% INPUTS:
%   calibration globals
%   settings -- column vectors in device settings

% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow gamma table dimensions to exceed device settings.
% 11/11/93   dhb   Update for new calData routines.
% 8/4/96     dhb   Update for stuff bag routines.
% 8/21/97    dhb   Update for structure.
% 4/5/02     dhb, ly  New calling interface.
% 8/3/07     dhb   Fix for [0-1] world.
% 5/08/14    npc   Modifications for accessing calibration data using a @CalStruct object.
%                  The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                  Passing a @CalStruct object is the preferred way because it results in 
%                  (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                  (b) better control over how the calibration data are accessed.

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Get gamma table
gammaTable = calStructOBJ.get('gammaTable');
gammaInput = calStructOBJ.get('gammaInput');
if (isempty(gammaTable))
	error('No gamma table present in calibration structure');
end

% Check dimensions and table sizes
[m,n] = size(settings);
[mg,ng] = size(gammaTable);
if (m > ng)
  error('Mismatch between primary coordinate dimensions and gamma table');
end

% Use a search routine to find the best gamma function
primary = zeros(m,n);
for i = 1:m
  [primary(i,:)] = SearchGammaTable(settings(i,:),gammaTable(:,i),gammaInput);
end

return

% This is the old OS 9 code, which doesn't work in the 0-1 world.

% Convert settings scale from [0:max-1] to [1:max]
%settings = settings+ones(m,n);

% Invert the gamma correction:  the settings happen to be the 
% indices into the gamma table.  Inverting the gamma correction 
% amounts to returning the elements of the gamma table corresponding
% to the settings.
%primary = zeros(m,n);
%for i = 1:m
%  primary(i,:) = gammaTable(settings(i,:),i)';
%end
