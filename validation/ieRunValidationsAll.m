function ieRunTutorialsAll
%ieRunTutorialsAll
%
% Description:
%   Run all of the isetbio validations that we think should work, and print out a report at the end
%   as to whether they threw errors, or not.  This does not check the validation data, it just
%   runs through the scripts.

% 10/19/17  dhb  Wrote this from tutorials script

% User/project specific preferences
p = struct(...
    'rootDirectory',            fileparts(which(mfilename())), ...
    'tutorialsSourceDir',       fullfile(isetbioRootPath, 'validation', 'scripts') ...                % local directory where tutorial scripts are located
    );

%% List of scripts to be skipped from automatic publishing.
%
% Anything with this in its path name is skipped.
scriptsToSkip = {...
    'codedevscripts' ...
    'xNeedChecking' ...
    'rgc' ...
    };

%% Use UnitTestToolbox method to do this.
UnitTest.runProjectTutorials(p, scriptsToSkip, 'All');
end