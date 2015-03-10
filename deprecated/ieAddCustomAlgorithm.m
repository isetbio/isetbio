function obj = ieAddCustomAlgorithm(obj,algType,handles,varargin)
%Deprecated - Add a custom computational algorithm to a pulldown list 
%
%   obj = ieAddCustomAlgorithm(obj,algType,handles,varargin)
%
% Users can add or delete a custom algorithms from several pulldown menus
% in ISET.  The selected algorithm is added to the vcSESSION file and from
% the GUI pull down list. Algorithms that can be modified in this way are
% currently stored in the vcimage window, oi window, and sensor window.
%
% The list of algorithms that are inserted and deleted is
%
%   Processor window:
%      {'demosaic','colordemosaic'}
%      {'illuminant correction'}
%      {'sensor conversion'}
%      {'processingmethod','render'}
%
%    OI window:
%      {'oicompute','oicomputemethod'}
%
%    Sensor window:
%      {'sensorcompute','sensorcomputemethod'} {'edgealgorithm'}
%
% To delete algorithms from the popup menus use ieDeleteCustomAlgorithm
%
% Example:
%    sensor = vcGetObject('isa');
%    handles = ieSessionGet('sensorhandles');
%    sensor = ieAddCustomAlgorithm(sensor,'sensorcompute',handles);
%    vcReplace ...
%
% Copyright ImagEval Consultants, LLC, 2005.

error('Deprecated');

switch lower(algType)
    
    case {'demosaic','colordemosaic'}
        % Obsolete
        currentAlg = get(handles.popDemosaic,'String'); 
    
        algFile = vcSelectAlgFile(['Processing',filesep,'Demosaic']);
        if isempty(algFile), return; end
        
        [p,baseName] = fileparts(algFile);
        if ismember(baseName,currentAlg)
            error('A function with this name is already in the list.');
        end
        obj = imageSet(obj,'demosaicmethod',baseName);
        
        newAlg = {baseName};
        ieSessionSet('adddemosaic',newAlg);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentAlg,newAlg);
        set(handles.popDemosaic,'String',newList);
        
    case {'colorbalance','balance'}
                % Obsolete

        currentAlg = get(handles.popBalance,'String'); 
    
        algFile = vcSelectAlgFile(['Processing',filesep,'Balance']);
        if isempty(algFile), return; end
        
        [p,baseName] = fileparts(algFile);
        if ismember(baseName,currentAlg)
            errordlg('A function with this name is already in the list.');
        end
        obj = imageSet(obj,'colorbalancemethod',baseName);
        
        newAlg = {baseName};
        ieSessionSet('addcolorbalance',newAlg);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentAlg,newAlg);
        set(handles.popBalance,'String',newList);    
        
    case {'conversion','colorconversion','colorconversionmethod'}
                % Obsolete

        currentAlg = get(handles.popColorConversionM,'String'); 
    
        algFile = vcSelectAlgFile(['Processing',filesep,'Conversion']);
        if isempty(algFile), return; end
        
        [p,baseName] = fileparts(algFile);
        if ismember(baseName,currentAlg)
            errordlg('A function with this name is already in the list.');
        end
        obj = imageSet(obj,'conversionmethod',baseName);
        
        newAlg = {baseName};
        ieSessionSet('addcolorconversion',newAlg);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentAlg,newAlg);
        set(handles.popColorConversionM,'String',newList);
        
    case {'processingmethod','render'}
                
        currentRender = get(handles.popCustomPath,'String'); 
    
        renderFile = vcSelectAlgFile(['Processing',filesep,'Render']);
        if isempty(renderFile), return; end
        
        [p,baseName] = fileparts(renderFile);
        if ismember(baseName,currentRender)
            errordlg('A function with this name is already in the list.');
        end
        obj = imageSet(obj,'rendermethod',baseName);
        
        renderFile = {baseName};
        ieSessionSet('addrendermethod',renderFile);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentRender,renderFile);
        set(handles.popCustomPath,'String',newList);
        
    case {'oicompute','oicomputemethod'}
        
        currentOI = get(handles.popCustom,'String'); 
        
        oiComputeFile = vcSelectAlgFile('OpticalImage');
        if isempty(oiComputeFile), return; end
        
        [p,baseName] = fileparts(oiComputeFile);
        if ismember(baseName,currentOI)
            set(handles.txtMessage,'string','An OI function with this name is already in the list.');
            return;
        end
        obj = oiSet(obj,'oiComputeMethod',baseName);
        
        oiComputeFile = {baseName};
        ieSessionSet('addOIComputemethod',oiComputeFile);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentOI,oiComputeFile);
        set(handles.popCustom,'String',newList);
        set(handles.txtMessage,'string','');

    case {'sensorcompute','sensorcomputemethod'}
        
        currentSensor = get(handles.popCustomCompute,'String'); 
    
        sensorComputeFile = vcSelectAlgFile('Sensor');
        if isempty(sensorComputeFile), return; end
        
        [p,baseName] = fileparts(sensorComputeFile);
        if ismember(baseName,currentSensor)
            errordlg('A function with this name is already in the list.');
        end
        obj = sensorSet(obj,'sensorComputeMethod',baseName);
        
        sensorComputeFile = {baseName};
        ieSessionSet('addsensorcompute',sensorComputeFile);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentSensor,sensorComputeFile);
        set(handles.popCustomCompute,'String',newList);
        
    case {'edgealgorithm'}
        
        currentEdge = get(handles.popCustomEdge,'String'); 
        
        edgeComputeFile = vcSelectAlgFile('Edge');
        if isempty(edgeComputeFile), return; end
        
        [p,baseName] = fileparts(edgeComputeFile);
        if ismember(baseName,currentEdge)
            errordlg('A function with this name is already in the list.');
        end
        obj = edgeSet(obj,'edgeoperator',baseName);
        
        edgeComputeFile = {baseName};
        ieSessionSet('addedgealgorithm',edgeComputeFile);
        
        % Set the popup list and select the new algorithm
        newList = cellMerge(currentEdge,edgeComputeFile);
        set(handles.popCustomEdge,'String',newList);
        
    otherwise
        error('Unknown algorithm type');
end

return;

%--------------------------------------
function fname = vcSelectAlgFile(algType)
%
%  fname = vcSelectDataFile(algType)
%
%Select a custom algorithm file.  This routine is really obsolete in a way.
%But I didn't have the energy to change all the stuff above.  So the
%routine is not really more like verifyAlgFile()

fname = [];

fullName = vcSelectDataFile('algorithm','r');
if isempty(fullName), 
    return;
else
    [p,fname,e] = fileparts(fullName);
    if ~strcmp(e,'.m'), errordlg('You must choose an m-file'); fname = []; end
    if isempty(which(fname)), warndlg('Cannot find file on the path'); end
end

return;
