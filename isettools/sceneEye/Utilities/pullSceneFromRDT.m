function success = pullSceneFromRDT(sceneName, scenePath)
% If a file isn't in the data folder, retrieve from RemoteDataToolbox
%
% Syntax:
%   success = pullSceneFromRDT(sceneName, scenePath)
%
% Description:
%    Check if given file exists in the data folder of iset3d. If it doesn't
%    try to pull it from the RemoteDataToolbox. The sceneName must match
%    the name of the artifact on RDT.
%
% Inputs:
%    sceneName - String. The artifact name to pull.
%    scenePath - String. Where to check if the file already exists .
%
% Outputs:
%    success   - Boolean. Whether or not the file (now) exists.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/17  TL   SCIEN Stanford, 2017
%    05/29/19  JNM  Documentation pass

% Download from RDT
if ~exist(scenePath, 'file')
    piPBRTFetch(sceneName, 'deletezip', true, 'pbrtversion', 3, ...
        'destination folder', fullfile(piRootPath, 'local', 'scenes'));
    % Check if file exists
    if ~exist(scenePath, 'file')
        success = false;
        error('Something went wrong when downloading the scene.')
    else
        % Success!
        fprintf('PBRT scene downloaded! File located at: %s\n', scenePath);
        success = true;
    end
else
    fprintf('Scene already exists in data folder. Skipping download.\n');
    success = true;
end

end
