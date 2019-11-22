function [oiObjects, seObjects] = downloadFromCloud(gcp, varargin)
% Download rendered data from the cloud. 
%
% Syntax:
%   [oiObjects, seObjects] = downloadFromCloud(gcp, [varargin])
%
% Description:
%    An variation of the function "downloadPBRT" from isetcloud but
%    adapated to work specifically for sceneEye. (We need to do some extra
%    processing specific to the sceneEye class.)
%
% Inputs:
%    gcp              - Object. The intalized gCloud object from isetcloud
%
% Outputs:
%    oiObjects        - Object. All of the optical images rendered.
%    seObjects        - Object. The corresponding scene eye objects. 
%
% Optional key/value pairs:
%    scaleIlluminance -  Boolean. If true, we scale the mean illuminance by
%                        the pupil diameter in piDat2ISET
%

% History:
%    04/26/18  TL   Created
%    05/29/19  JNM  Documentation pass

%% Initialize
p = inputParser;
p.addRequired('gcp', @(x)(isa(x, 'gCloud')));
p.addParameter('scaleIlluminance', true, @islogical);

p.parse(gcp, varargin{:});
scaleIlluminance = p.Results.scaleIlluminance;

oiObjects = [];
seObjects = [];

%% Download
oiObjects = gcp.downloadPBRT(gcp.miscDescriptor(1).recipe, ...
    'scaleIlluminance', scaleIlluminance);

for ii = 1:length(oiObjects)
    % Get corresponding sceneEye object
    % Check, is this correct?
    seObjects{ii} = gcp.miscDescriptor(ii);

    % Set the parameters correctly for the optical image
    if seObjects{ii}.debugMode == 1
        % seObject is a scene. Don't try to set it as an optical image,
        % just return it. 
        oiObjects{ii} = oiObjects{ii};
    else
        oiObjects{ii} = seObjects{ii}.setOI(oiObjects{ii}, ...
            'scaleIlluminance', scaleIlluminance);
    end
end

end
