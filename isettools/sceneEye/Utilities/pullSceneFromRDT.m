function success = pullSceneFromRDT(sceneName,scenePath)
%DOWNLOADSCENEFROMRDT Check if given file exists in the data folder of
%   iset3d. If it doesn't try to pull it from the RemoteDataToolbox.
% sceneName must match the name of the artifact on RDT. 
%
% Input
%   sceneName - artifact name to pull
%   scenePath - where to check if the file already exists 
%
% TL, SCIEN Stanford, 2017

% Download from RDT
if(~exist(scenePath,'file'))
    piPBRTFetch(sceneName,'deletezip',true,...
        'pbrtversion',3,...
        'destination folder',fullfile(piRootPath,'data','V3'));
    % Check if file exists
    if(~exist(scenePath,'file'))
        error('Something went wrong when downloading the scene.')
        success = false;
    else
        % Success!
        fprintf('PBRT scene downloaded! File is located at: %s \n',scenePath);
        success = true;
    end
    
else
    fprintf('Scene already exists in data folder. Skipping download.\n');
    success = true;
end

            

end

