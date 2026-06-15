function testResults = isetbioExamplesTest()
% isetbioExamplesTest - Runs all s_* and t_* scripts in the examples directory securely.
%
% Usage:
%   testResults = isetbioExamplesTest;
%
% It evaluates each script in the base workspace so that `clear` commands
% within the scripts do not destroy the state of this runner function.
% It catches errors and prints a summary report at the end.

    % Find the Examples directory
    examplesDir = fullfile(isetRootPath, 'examples');
    
    % Recursively find all scripts starting with s_ or t_
    % Wait, this is ISETBio, so root is isetbioRootPath, not isetRootPath
    % Let's fix that!
    
