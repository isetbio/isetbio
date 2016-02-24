% Utility to remove one validation ground truth data set (both fast and full)
%
% Usage: ieDeleteValidationFile('osBiophysObject')
%

function ieCompareLocalToRemoteValidationFile(validationFile)

    list = rdtListLocalArtifacts(...
        getpref('isetbioValidation', 'remoteDataToolboxConfig'), ...
        'validation/full', 'artifactId', validationFile, 'type', 'mat');

    if isempty(list)
        fprintf('Did not find any artifacts matching ''%s''.', validationFile)
    else
        fprintf('Found the following local/remote artifacts: \n');
        for k = 1:numel(list)
            fprintf('[%02d].  ''%s''  remotePath:%s   localPath: %s  \n', k, list(k).artifactId, list(k).remotePath, list(k).localPath);
        end
        
        fprintf('Local file contents\n')
        load(list(1).localPath)
        whos
        validationData
        hostInfo
%         validationData.diffractionOI
%         validationData.diffractionOI.data
%         validationData.diffractionOI.optics
%         validationData.diffractionOI.diffuser
        
        whos
        pause
        client = RdtClient('isetbio');
        client.crp(list(k).remotePath);
        %client.openBrowser();
        [artifactData, artifactInfo] = client.readArtifact(validationFile, 'type', 'mat');
        artifactData.validationData
        artifactData.hostInfo
        
        clear all
        load('/Users/nicolas/Downloads/test.mat')
        whos
        validationData
        hostInfo
%         validationData.diffractionOI
%         validationData.diffractionOI.data
%         validationData.diffractionOI.optics
%         validationData.diffractionOI.diffuser
    end
end