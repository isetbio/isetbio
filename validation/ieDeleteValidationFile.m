function ieDeleteValidationFile
% ieDeleteValidationFile
%
%% Utility to remove one validation ground truth data set (both fast and full)

% Which project
thisProject = 'isetbio';
validationFileToBeDeleted =  selectValidationFile(thisProject);

list = rdtListLocalArtifacts(...
    getpref(thisProject, 'remoteDataToolboxConfig'), ...
    'validation', 'artifactId', validationFileToBeDeleted);

if isempty(list)
    fprintf('Did not find any artifacts matching ''%s''.\n', validationFileToBeDeleted)
else
    fprintf('Found the following local/remote artifacts: \n');
    for k = 1:numel(list)
        fprintf('[%02d].  ''%s''  remotePath:%s   localPath: %s  \n', k, list(k).artifactId, list(k).remotePath, list(k).localPath);
    end
    fprintf('Hit enter to delete these %d artifacts. ', numel(list))
    pause;
    client = RdtClient(getpref(thisProject,'remoteDataToolboxConfig'));
    
    % Log into archiva
    if isempty(client.configuration.password) ...
            && ~isempty(client.configuration.username)
        client.credentialsDialog();
    end
    
    for k = 1:numel(list)
        fprintf('Deleting %s/%s\n', list(k).remotePath, list(k).artifactId);
        deleted = rdtDeleteLocalArtifacts(client.configuration, list(k));
        deleted = rdtDeleteArtifacts(client.configuration, list(k));
    end
end
end

function validationFileToBeDeleted = selectValidationFile(thisProject)

UnitTest.usePreferencesForProject(thisProject);

validationFileToBeDeleted = UnitTest.selectScriptFromExistingOnes('prompt','Enter script no. to delete data for that script');
idx = strfind(validationFileToBeDeleted,'/');
validationFileToBeDeleted = validationFileToBeDeleted(idx(end)+1:end);
validationFileToBeDeleted = strrep(validationFileToBeDeleted, 'v_', '');
validationFileToBeDeleted = strrep(validationFileToBeDeleted, '.m', '');
end