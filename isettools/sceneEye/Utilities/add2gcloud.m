function [cloudFolder,zipFileName] = add2gcloud(gcp,se,varargin)
%ADD2GCLOUD Add a sceneEye object to a gcp object so that it can be
%rendered on gCloud.
%
% isetcloud currently works specifically with iset3d recipes. The sceneEye
% class, however, is different than the recipe class. (From the point of
% view of sceneEye, the recipe is the "raw" data and functions within
% sceneEye work to change parts of the recipe.)
%
% This function does all the necessary steps to add a sceneEye's internal
% recipe to a gCloud object.
%
% Inputs:
%    gcp - the intalized gCloud object from isetcloud
%    sceneEye - the sceneEye class to add to the gcp targets.
%    uploadZip - whether the folder should be zipped up and uploaded to
%               gCloud at this point. Users shoudl call add2gcloud with
%               this flag when the final sceneEye class is being uploaded.
%    varargin  - An optional length of key/value pairs describing the scene
%
% Outputs:
%    success - success flag?
%
% History:
%    4/23/18  TL   Created
%%
p = inputParser;
p.addRequired('gcp',@(x)(isa(x,'gCloud')));
p.addRequired('se',@(x)(isa(x,'sceneEye')));
p.addParameter('uploadZip',false,@islogical);

p.parse(gcp,se,varargin{:});

uploadZip = p.Results.uploadZip;

% Empty until the zip is uploaded.
cloudFolder = [];
zipFileName = [];

%%

% Write the sceneEye recipe out
seNew = se.write();

% Upload the recipe 
if(uploadZip)
    [cloudFolder, zipFileName] = gcp.uploadPBRT(seNew.recipe,'overwrite zip',true);
else
    gcp.uploadPBRT(seNew.recipe,'upload zip',false,'overwrite zip',false);
end

% Add the new target operation
addPBRTTarget(gcp,seNew.recipe);


end

