% The daylight files in isettools/data/lights
%
% Description:
%    A script covering the use of the daylight files in
%    isettools/data/lights.
%

% History:
%    XX/XX/XX  BW   Created by Wandell
%    11/23/18  JNM  Formatting

%%
chdir(fullfile(isetRootPath, 'data', 'lights'))

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

%% daylightPsychBldg
% A single daylight example measured at Stanford Psychology

% I know we have a large data set of 10,000 measurements.  I am not
% sure why we have this single example here.
[daylightExamples, wave, comment] = ieReadSpectra('DaylightPsychBldg');
disp(comment)

vcNewGraphWin;
plot(wave, daylightExamples);
grid on;
xlabel('Wave');
ylabel('Relative energy');