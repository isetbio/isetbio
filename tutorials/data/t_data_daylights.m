% The daylight files in ISETCam's data/lights directory
%
% Description:
%    A script covering the use of the daylight files in
%    ISETCam's data/lights directory.
%

% History:
%    XX/XX/XX  BW   Created by Wandell
%    11/23/18  JNM  Formatting

%% The data live in ISETCam's data/lights/daylight directory
% ieReadSpectra finds them on the MATLAB path, so there is no need to
% change the working directory to read them.

%% cieDaylightBasis
% Nov 2017 BW removed cieDaylight and cieDay which were lower resolution
% representations of the same data
% I think the original source for these are the PsychTOolbox.

wave = 380:1:780;
[cieDay, ~, comment] = ieReadSpectra('cieDaylightBasis', wave);
disp(comment)

vcNewGraphWin;
plot(wave, cieDay);
grid on;
xlabel('Wave');
ylabel('Relative energy');

%% daylightStanford
% Daylight measured at Stanford by J. DiCarlo, winter 2000. This is the
% large measured set, about a thousand spectra, that the single
% 'DaylightPsychBldg' example used to stand in for. That file was removed
% when ISETCam reorganized data/lights, so read the full set instead.
[daylightExamples, wave, comment] = ieReadSpectra('daylightStanford');
disp(comment)

% Plot a readable subsample rather than all of the measurements.
vcNewGraphWin;
plot(wave, daylightExamples(:, 1:25:end));
grid on;
xlabel('Wave');
ylabel('Relative energy');
title(sprintf('%d of %d daylight measurements', ...
    numel(1:25:size(daylightExamples, 2)), size(daylightExamples, 2)));
