% Utility to compare local to remote validationd data
%
% Usage: ieCompareLocalToRemoteValidationFile('oi')
%

function ieCompareLocalToRemoteValidationFile(validationFile)

    list = rdtListLocalArtifacts(...
        getpref('isetbio', 'remoteDataToolboxConfig'), ...
        'validation/full', 'artifactId', validationFile, 'type', 'mat');

    if isempty(list)
        fprintf('Did not find any artifacts matching ''%s''.', validationFile)
    else
        fprintf('Found the following local/remote artifacts: \n');
        for k = 1:numel(list)
            fprintf('[%02d].  ''%s''  remotePath:%s   localPath: %s  \n', k, list(k).artifactId, list(k).remotePath, list(k).localPath);
        end
        
        fprintf('Local file contents (file: %s)\n', list(1).localPath)
        load(list(1).localPath)
        validationData
        hostInfo
        fprintf('Hit enter to fetch remote data \n');
        pause
        
        client = RdtClient('isetbio');
        client.crp(list(k).remotePath);
        client.openBrowser();
        [artifactData, artifactInfo] = client.readArtifact(validationFile, 'type', 'mat');
        artifactData.validationData
        artifactData.hostInfo
        pause
        
      
    end
end