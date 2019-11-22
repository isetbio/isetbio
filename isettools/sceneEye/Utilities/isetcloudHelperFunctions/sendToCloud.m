function [cloudFolder, zipFileName] = sendToCloud(gcp, se, varargin)
% Add a sceneEye object to gcp object so that it can be rendered on gCloud
%
% Syntax:
%   [cloudFolder, zipFileName] = sendToCloud(gcp, se, [varargin])
%
% Description:
%    isetcloud currently works specifically with iset3d recipes. The
%    sceneEye class, however, is different than the recipe class. (From the
%    point of view of sceneEye, the recipe is the "raw" data and functions
%    within sceneEye work to change parts of the recipe.)
%
%    This function does all the necessary steps to add sceneEye object's
%    internal recipe to a gCloud object.
%
% Inputs:
%    gcp         - Object. The intalized gCloud object from isetcloud
%    sceneEye    - Object. The sceneEye class to add to the gcp targets.
%
% Outputs:
%    cloudFolder - String. The location where things are sent on gcloud.
%    zipFileName - String. The filename for the resource zipfile to be sent
%                  to the cloud.
%
% Optional key/value pairs:
%    uploadZip   - Boolean. Whether the folder should be zipped up and
%                  uploaded to gCloud at this point. Users should call
%                  sendToCloud with this flag only when they are uploading
%                  their final sceneEye object.
%

% History:
%    04/23/18  TL   Created
%    05/29/19  JNM  Documentation pass

%% Initialize
p = inputParser;
p.addRequired('gcp', @(x)(isa(x, 'gCloud')));
p.addRequired('se', @(x)(isa(x, 'sceneEye')));
p.addParameter('uploadZip', false, @islogical);

p.parse(gcp, se, varargin{:});
uploadZip = p.Results.uploadZip;

% Empty until the zip is uploaded.
cloudFolder = [];
zipFileName = [];

%%
% Write the sceneEye recipe out
seNew = se.write();

% Upload the recipe 
if(uploadZip)
    [cloudFolder, zipFileName] = ...
        gcp.uploadPBRT(seNew.recipe, 'overwrite zip', true);
else
    gcp.uploadPBRT(seNew.recipe, 'materials', false, ...
        'geometry', false, 'resources', false, 'overwrite zip', false);
end

% Add the new target operation
addPBRTTarget(gcp, seNew.recipe);

% We're going to keep track of the sceneEye objects for each target by
% adding another variable to gcp.
gcp.miscDescriptor = cat(1, gcp.miscDescriptor, seNew.copy);

end
