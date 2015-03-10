function obj = ieDeleteCustomAlgorithm(obj,algType,handles)
% Deprecated - Gateway routine for deleting a custom algorithm from a pulldown list.
%
%   obj = ieDeleteCustomAlgorithm(obj,algType,handles)
%
%Purpose:
%  Users can delete a custom algorithms from several pulldown menus in
%  ISET.  The selected algorithm is removed from the vcSESSION file and
%  removed from the GUI pull down list. Algorithms that can be modified in
%  this way are currently stored in the vcimage window, oi window, and
%  sensor window.  The list is
%
%   Processor window:
%      {'demosaic','colordemosaic'}
%      {'colorbalance','balance'}
%      {'conversion','colorconversion','colorconversionmethod'}
%      {'processingmethod','render'}
%
%    OI window:
%      {'oicompute','oicomputemethod'}
%
%    Sensor window:
%      {'sensorcompute','sensorcomputemethod'} {'edgealgorithm'}
%
% Algorithms are added using ieAddCustomAlgorithm
%
% Example:
%    vci = ieDeleteCustomAlgorithm(vci,algType,handles)
%
% Copyright ImagEval Consultants, LLC, 2003.

error('Deprecated.')

if notDefined('obj'), error('Object must be defined.'); end
if notDefined('algType'), error('Algorithm type must be defined'); end

pauseTime = 2;  % sec

% If the handles are not passed in, we can get them this way.
if notDefined('handles')
    if ~checkfields(obj,'type')
        error('Object type is not defined.'); 
    else
        objType = vcEquivalentObjtype(obj.type);
        switch lower(objType)
            case 'opticalimage'
                handles = ieSessionGet('opticalimagehandle');
            case 'isa'
                handles = ieSessionGet('sensorimagehandle');
            case 'vcimage'
                handles = ieSessionGet('vcimagehandle');
            otherwise
                error('Object type not managed by this routine.');
        end
    end
end
            
switch lower(algType)
    
    case {'demosaic','colordemosaic'}
        defaultAlg = vcAlgorithms('demosaic');
        
        customAlg = ieSessionGet('demosaiclist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Demosaics');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('demosaiclist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = imageSet(obj,'demosaicmethod',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popDemosaic,'String',newList);
        
    case {'colorbalance','balance'}
        defaultAlg = vcAlgorithms('balance');
        
        customAlg = ieSessionGet('colorbalancelist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Balance');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('colorbalancelist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = imageSet(obj,'colorbalancemethod',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popBalance,'String',newList);
        
    case {'colorconversion','conversion','colorconversionmethod'}
        defaultAlg = vcAlgorithms('colorconversion');
        
        customAlg = ieSessionGet('colorconversionlist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('colorconversionlist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = imageSet(obj,'colorconversionmethod',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popColorConversionM,'String',newList);
        
    case {'render'}
        defaultAlg = vcAlgorithms('render');
        
        customAlg = ieSessionGet('renderlist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('renderlist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = imageSet(obj,'render',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popCustomPath,'String',newList);
        
    case {'oicompute','oicomputemethod'}
        % From OI window ... delete some of the custom oiCompute options
        defaultAlg = vcAlgorithms('oicompute');
        
        customAlg = ieSessionGet('oicomputelist');
        if isempty(customAlg)
            set(handles.txtMessage,'string','No custom algorithms.',handles,pauseTime); 
            return; 
        end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('oicomputelist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = oiSet(obj,'oiCompute',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popCustom,'String',newList);
        set(handles.txtMessage,'string','');

    case {'sensorcompute','sensorcomputemethod'}
        % From Sensor window ... delete some of the custom sensorCompute options
        defaultAlg = vcAlgorithms('sensorcompute');
        
        customAlg = ieSessionGet('sensorcomputelist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('sensorcomputelist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = sensorSet(obj,'sensorCompute',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popCustomCompute,'String',newList);
        
    case {'render'}
        defaultAlg = vcAlgorithms('render');
        
        customAlg = ieSessionGet('renderlist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('renderlist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = imageSet(obj,'render',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popCustomPath,'String',newList);
        
    case {'edgealgorithm'}
        % From Sensor window ... delete some of the sensorCompute options
        defaultAlg = vcAlgorithms('edgealgorithms');
        
        customAlg = ieSessionGet('edgealgorithmlist');
        if isempty(customAlg), ieInWindowMessage('No custom algorithms.',handles,pauseTime); return; end
        
        [dList,OK] = listdlg('ListString',customAlg,'Name','Delete Conversion');
        if ~OK, return; end
        
        customAlg = cellDelete(customAlg,dList);
        ieSessionSet('edgealgorithmlist',customAlg);
        
        % Set the default to the first on the list, I suppose.  Or find
        % another process.
        obj = edgeSet(obj,'edgealgorithm',defaultAlg{1});
        
        newList = cellMerge(defaultAlg,customAlg);
        set(handles.popCustomEdge,'String',newList);
        
    otherwise
        error('Unknown algorithm type');
end
    
end
