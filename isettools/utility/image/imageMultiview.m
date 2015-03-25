function selectedObjs = imageMultiview(objType, selectedObjs, singlewindow)
% Display multiple images of selected GUI objects
%
%  selectedObjs = imageMultiview(objType, whichObjs, singlewindow)
%
% This routine lets the user compare the images side by side, rather than
% flipping through them in the GUI window.
%
% objType:       Which window (scene, oi, sensor, or vcimage)
% selectedObjs:  List of the selected object numbers, e.g., [1 3 5]
% singlewindow:  Put the images in subplots of a single figure (true) or in
%                different figures (default = false);
%
% See also: imageMontage
%
% Example:
%  objType = 'scene';
%  imageMultiview(objType);
%
%  selectedObjs = [1 6];
%  imageMultiview(objType,whichObj);
%
%  objType = 'vcimage';
%  selectedObjs = [2 3 5];
%  imageMultiview(objType,whichObj, true);
%
% Copyright Imageval Consultants, LLC, 2013

if notDefined('objType'), error('Object type required.'); end
if notDefined('singlewindow'), singlewindow = false; end

% Allows some aliases to be used
objType = vcEquivalentObjtype(objType);

% Get the objects
[objList, nObj] = vcGetObjects(objType);
if  isempty(objList)
    fprintf('No objects of type %s\n',objType);
    return;
end

% Show a subset or all
if notDefined('selectedObjs')
    lst = cell(1,nObj);
    for ii=1:nObj, lst{ii} = objList{ii}.name; end
    selectedObjs = listdlg('ListString',lst);
end

% Set up the subplots or multiple window conditions
if singlewindow
    if nObj > 3
        rWin = ceil(sqrt(length(selectedObjs)));
        cWin = rWin; fType = 'upper left';
    else
        rWin = nObj; cWin = 1; fType = 'tall';
    end
else   rWin = []; fType = 'upper left';
end
gam = 1;  % Figure out a rationale for this.
subCount = 1; % Which subplot are we in

%% This is the display loop
for ii=selectedObjs
    if (~singlewindow || subCount == 1), f = vcNewGraphWin([],fType); end
    if singlewindow,   subplot(rWin,cWin,subCount); subCount = subCount+1; end
    switch objType
        case 'SCENE'
            sceneShowImage(objList{ii},true,gam);
            t = sprintf('Scene %d - %s',ii,sceneGet(objList{ii},'name'));
            
        case 'OPTICALIMAGE'
            oiShowImage(objList{ii},true,gam);
            t =sprintf('OI %d - %s',ii,oiGet(objList{ii},'name'));
            
        case 'ISA'
            scaleMax = 1;
            sensorShowImage(objList{ii},gam,scaleMax);
            t = sprintf('Sensor %d - %s',ii,sensorGet(objList{ii},'name'));
            
        case 'VCIMAGE'
            imageShowImage(objList{ii},gam,true,f);
            t = sprintf('VCI %d - %s',ii,imageGet(objList{ii},'name'));
            
        otherwise
            error('Unsupported object type %s\n', objType);
    end
    
    % Label the image or window
    if singlewindow,      title(t)
    else                  set(gcf,'name',t);
    end
    
end

end