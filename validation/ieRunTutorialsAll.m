function ieRunTutorialsAll
%ieRunTutorialsAll
%
% Syntax
%    ieRunTutorialsAll
%
% Description
%   Run all of the isetbio tutorials that we think should work, and print
%   out a report at the end as to whether they threw errors, or not.
%   Scripts inside of isetbioRootPath/tutorials are run, except that
%   scripts within directories 'underDevelopment' or 'support' in
%   their names are skipped.
%
% 
% 07/26/17  dhb  Wrote this, because we care.

% User/project specific preferences
p = struct(...
    'rootDirectory',            fileparts(which(mfilename())), ...
    'tutorialsSourceDir',       fullfile(isetbioRootPath, 'tutorials') ...                % local directory where tutorial scripts are located
    );

%% List of scripts to be skipped from automatic publishing.
%
% Anything with this in its path name is skipped.
scriptsToSkip = {...
    'underDevelopment' ...
    'support' ...
    't_peripheralOpticsAndMosaic' ...'
    't_cMosaicBenchMark.m' ...
    };

%% Use UnitTestToolbox method to do this.
UnitTest.runProjectTutorials(p, scriptsToSkip, 'All');
end