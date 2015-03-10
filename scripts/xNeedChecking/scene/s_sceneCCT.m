%% s_sceneCCT
%
% Calculate the correlated color temperature of a scene illuminant
%
% The correlated color temperature is an assessment of the yellow-blue
% appearance of the illuminant.  Low temperatures (3500K) are yellowish and
% high temperatures (7000K)are bluish.  A typical mix of blue sky and the
% sun is 5500 or 6500K.
%
% The CCT calculation was developed by Wyszecki and Judd in the 1960s.  It
% relies on the old (u,v) format, not the (u',v') format.  See notes in
% xyz2uv.
%
% Copyright Imageval Consulting, LLC 2013

%%
s_initISET

%% Create an energy representation of the illuminant

wave = 400:5:720;
spd = blackbody(wave,3500);

vcNewGraphWin; 
plot(wave,spd);
grid on; xlabel('Wavelength (nm)'); ylabel('Energy (watts/sr/nm/m^2)')

%% spd2cct converts energy to XYZ and then uv (not uprime,vprime)

cTemp = 3500;
d = blackbody(wave, cTemp);
fprintf('Estimated CCT %.1f and actual %.1f\n',spd2cct(wave,d),cTemp)

cTemp = 6500;
d = blackbody(wave, cTemp);
fprintf('Estimated CCT %.1f and actual %.1f\n',spd2cct(wave,d),cTemp)

cTemp = 8500;
d = blackbody(wave, cTemp);
fprintf('Estimated CCT %.1f and actual %.1f\n',spd2cct(wave,d),cTemp)

%% Or, do several spds at once

cTemps = 4500:1000:8500;
spd = blackbody(wave, cTemps);

vcNewGraphWin; 
plot(wave,spd);
grid on; xlabel('Wavelength (nm)'); ylabel('Energy (watts/sr/nm/m^2)')

%% Print out the list
disp('-----')
for ii=1:length(cTemps)
    fprintf('Estimated CCT %.1f and actual %.1f\n',spd2cct(wave,spd(:,ii)),cTemps(ii));
end

%% Notes
%
% This is the routine we use to name an illuminant in the scene, when there
% is no name given by the user.
%

%% END

