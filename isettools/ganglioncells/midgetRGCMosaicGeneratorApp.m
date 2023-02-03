classdef midgetRGCMosaicGeneratorApp < handle


    properties  (GetAccess=public, SetAccess=private)
        % GUI components
        mainView;

        % Current action to perform
        currentAction;

        % State (mosaic and optics params)
        simulation;
    end

    properties (GetAccess=private, SetAccess=private)
        theMosaicGeometryTable;
        theOpticsTable;
        theRTVFmodelTable;
    end

    methods
        % Constructor
        function obj = midgetRGCMosaicGeneratorApp()
            
            generateGUI(obj);
            initializeState(obj);

        end
    end
end


function initializeState(obj)

     theParamsFile = midgetRGCMosaicGenerator.paramsFile();

     if (isfile(theParamsFile))
         load(theParamsFile, ...
             'mosaicCenterParams', ...
             'rfModelParams', ...
             'opticsParams');
         obj.simulation.mosaicCenterParams = mosaicCenterParams;
         obj.simulation.rfModelParams  = rfModelParams;
         obj.simulation.opticsParams = opticsParams;
     else
        % Generate the default params structs
        [obj.simulation.mosaicCenterParams, ...
         obj.simulation.rfModelParams, ...
         obj.simulation.opticsParams] = midgetRGCMosaicGenerator.generateDefaultMosaicAndOpticsParamStructs();
     end

     % Update theMosaicGeometryTable
     d = {};
     d{1} = obj.simulation.mosaicCenterParams.positionDegs(1);
     d{2} = obj.simulation.mosaicCenterParams.positionDegs(2);
     d{3} = obj.simulation.mosaicCenterParams.sizeDegs(1);
     d{4} = obj.simulation.mosaicCenterParams.sizeDegs(2);
     set(obj.theMosaicGeometryTable,'Data',d)

     d = {};
     d{1} = obj.simulation.opticsParams.ZernikeDataBase;
     d{2} = obj.simulation.opticsParams.subjectRankOrder;
     d{3} = obj.simulation.opticsParams.pupilDiameterMM;
     set(obj.theOpticsTable, 'Data', d);

     d = {};
     d{1} = obj.simulation.rfModelParams.H1cellIndex;
     d{2} = obj.simulation.rfModelParams.eccentricitySamplingGridHalfSamplesNum;
     set(obj.theRTVFmodelTable, 'Data', d)
     drawnow
end


function executeButtonAction(btn, app)

    switch app.currentAction
        case "compute: center-connected mRGC mosaic"
            midgetRGCMosaicGenerator.generateCenterConnectedMosaic(...
                app.simulation.mosaicCenterParams);

        case "compute: all R2VFT objects"

            RTVobjIndicesToBeComputed = 'all';

            % If we had a crash, compute remaining RTVFobjects
            % Then update the list with the ones that were computed before
            % the crash, using 
            % "compute: manually replace a specific R2VFT object" for each
            % of these pre-crash computed RTVF objects
            %RTVobjIndicesToBeComputed = [22:22];
            
            RTVobjIndicesToBeComputed 
            pause
            
            midgetRGCMosaicGenerator.generateR2VFTobjects(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                'RTVobjIndicesToBeComputed', RTVobjIndicesToBeComputed);

        case "inspect: single R2VFT object file"
            midgetRGCMosaicInspector.quicklyInspectSingleRTVFobjectFile();

        case "inspect: all R2VFT objects file"
            midgetRGCMosaicInspector.quicklyInspectAllRTVFobjectsFile();

        case "compute: refit R2VFT object(s)"
            % Ask user which R2VFT objects to refit 
            % (single position,  multiple center cones num)
            targetPosition = [];
            while (numel(targetPosition) ~= 2)
                targetPosition = input('Enter position for which to update the RTVF object ([x y]): ');
            end
            targetRFcenterConesNum = input('Enter center cones num for which to update the RTVF object (e.g, 1, 2): ');
            targetRFcenterConeType = input('Enter center cone type for which to update the RTVF object (L or M or hit Enter for both): ', 's');

            % Re-generate R2VFT objects at a specific position
            midgetRGCMosaicGenerator.generateR2VFTobjects(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                'updateRTVFobjectAtPosition', targetPosition, ...
                'updateRTVFobjectWithCenterConesNum', targetRFcenterConesNum, ...
                'updateRTVFobjectWithCenterConeType', targetRFcenterConeType)

        case "compute: manually replace a specific R2VFT object"
            midgetRGCMosaicGenerator.replaceSpecificR2VFTobject()

        case "export: retinal cone pooling params for all fitted R2VFT objects"
            midgetRGCMosaicInspector.exportRetinalConePoolingParamsForAllFittedRTVFTobjects();

        case "compute: center-surround cone pooling kernels"
            midgetRGCMosaicGenerator.generateCenterSurroundConePoolingKernels(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams)

        case "visualize: spatial RFs"
            % Ask user how many RGCs to visualize
            maxRGCsNum = input('How many RGC RFs to visualize? ([Hit enter for all): ');

            midgetRGCMosaicInspector.visualizeSpatialRFmaps(...
                app.simulation.mosaicCenterParams, ...
                maxRGCsNum);
         
        case "compute: frozen midgetRGCMosaic (with current center-surround cone pooling weights)"
            midgetRGCMosaicGenerator.freezeMosaic(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams)

        case "validate: pre-compute LM-non-opponent cone mosaic STFs"
            % Ask user whether to use parfor or not
            useParfor = input('Use parfor for the cone mosaic computes? ([Hit enter to use parfor): ', 's');
            if (isempty(useParfor)) || (ischar(useParfor) && strcmpi(useParfor, 'y'))
                useParfor = true;
            else
                useParfor = false;
            end
            
            midgetRGCMosaicInspector.preComputeConeMosaicLMnonOpponentSTFs(...
                 app.simulation.mosaicCenterParams, ...
                 app.simulation.rfModelParams, ...
                 app.simulation.opticsParams, ...
                 useParfor);

        case "validate: compute LM-non-opponent RGC STFs (midgetRGCMosaic object handles the computation of cone mosaic responses)"
            % Ask user whether to use parfor or not
            useParfor = input('Use parfor for the midgetRGCMosaic computes? ([Hit enter to use parfor): ', 's');
            if (isempty(useParfor)) || (ischar(useParfor) && strcmpi(useParfor, 'y'))
                useParfor = true;
            else
                useParfor = false;
            end
            
            midgetRGCMosaicInspector.computeMosaicLMnonOpponentSTFs(...
                 app.simulation.mosaicCenterParams, ...
                 app.simulation.rfModelParams, ...
                 app.simulation.opticsParams, ...
                 useParfor);

        case "validate: compute LM-non-opponent RGC STFs from pre-computed LM-non-opponent cone mosaic STFs"
            useParfor = true;
            midgetRGCMosaicInspector.computeMosaicLMnonOpponentSTFsFromConeMosaicLMnonOpponentSTFs(...
                 app.simulation.mosaicCenterParams, ...
                 app.simulation.rfModelParams, ...
                 app.simulation.opticsParams, ...
                 useParfor);

        case "validate: fit LM-non-opponent STFs"
            % Ask user how many RGCs to analyze
            maxRGCsNum = input('How many RGCs to fit? ([Hit enter for all): ');

            midgetRGCMosaicInspector.fitMosaicSTFs(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                maxRGCsNum);


        case "validate: visualize STF fits"
            midgetRGCMosaicInspector.visualizeMosaicSTFfits(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams);

        otherwise
            error('Unknown action: ''%s''.', app.currentAction);
    end
end

function dropDownActionChanged(src,event, app)
    app.currentAction = event.Value;
end

function generateGUI(obj)

    % Create figure window
    obj.mainView = uifigure('Position', [30 1500 1000 360], ...
        'WindowStyle', 'AlwaysOnTop', ...
        'Scrollable', 'on', ...
        'Color', [0.3 0.3 0.3], ...
        'Resize', 'off');

    obj.mainView.Name = "midgetRGCMosaic generator & inspector";

    % Manage app layout
    layoutRows = 5;
    layoutCols = 2;
    gl = uigridlayout(obj.mainView,[layoutRows layoutCols]);
    gl.BackgroundColor = [0.3 0.3 0.3];
    gl.RowHeight = {'1x', '1x', '1x', '1x', '1x'};
    gl.ColumnWidth = {'fit','6x'};

    theActionLabel  = uilabel(gl);
    theMosaicGeometryButton  = uibutton(gl);
    obj.theMosaicGeometryTable = uitable(gl);
    obj.theOpticsTable = uitable(gl);
    obj.theRTVFmodelTable = uitable(gl);
    theOpticsButton = uibutton(gl);
    theRTVFmodelButton = uibutton(gl);
    theActionDropdown = uidropdown(gl);
    theExecuteActionButton = uibutton(gl);
    %theExitButton = uibutton(gl);

    % TheMosaicGeometry button
    theMosaicGeometryButton.Layout.Row = 3;
    theMosaicGeometryButton.Layout.Column = 1;
    theMosaicGeometryButton.Text = "mosaic params";
    theMosaicGeometryButton.FontSize = 20;
    theMosaicGeometryButton.FontWeight = 'Bold';
    theMosaicGeometryButton.BackgroundColor = [0.4 0.4 0.4];
    theMosaicGeometryButton.FontColor = [0.3 1 0.8];
    theMosaicGeometryButton.ButtonPushedFcn = @(btn,event) updateMosaicParams(btn, obj);

    obj.theMosaicGeometryTable.Layout.Row = 3;
    obj.theMosaicGeometryTable.Layout.Column = 2;
    obj.theMosaicGeometryTable.RowName = {};
    obj.theMosaicGeometryTable.ColumnName = {'x-position (degs)', 'y-position (degs)', 'x-size (degs)', 'y-size (degs)'};
    obj.theMosaicGeometryTable.Data = {};
    obj.theMosaicGeometryTable.FontSize = 18;
    obj.theMosaicGeometryTable.BackgroundColor = [0.4 0.4 0.4];
    addStyle(obj.theMosaicGeometryTable,uistyle('HorizontalAlignment','center', 'FontColor', [0.3 1 0.8]))

    % TheOptics button
    theOpticsButton.Layout.Row = 4;
    theOpticsButton.Layout.Column = 1;
    theOpticsButton.Text = "optics params";
    theOpticsButton.FontSize = 20;
    theOpticsButton.FontWeight = 'Bold';
    theOpticsButton.BackgroundColor = [0.4 0.4 0.4];
    theOpticsButton.FontColor = [0.3 1 0.8];
    theOpticsButton.ButtonPushedFcn = @(btn,event) updateOpticsParams(btn, obj);

    obj.theOpticsTable.Layout.Row = 4;
    obj.theOpticsTable.Layout.Column = 2;
    obj.theOpticsTable.RowName = {};
    obj.theOpticsTable.ColumnName = {'Zernike Database', 'subject rank', 'pupil diameter (mm)'};
    obj.theOpticsTable.Data = {};
    obj.theOpticsTable.FontSize = 18;
    obj.theOpticsTable.BackgroundColor = [0.4 0.4 0.4];
    addStyle(obj.theOpticsTable, uistyle('HorizontalAlignment','center', 'FontColor', [0.3 1 0.8]))


    % theRTVFmodel button
    theRTVFmodelButton.Layout.Row = 5;
    theRTVFmodelButton.Layout.Column = 1;
    theRTVFmodelButton.Text = "optics params";
    theRTVFmodelButton.FontSize = 20;
    theRTVFmodelButton.FontWeight = 'Bold';
    theRTVFmodelButton.BackgroundColor = [0.4 0.4 0.4];
    theRTVFmodelButton.FontColor = [0.3 1 0.8];
    theRTVFmodelButton.ButtonPushedFcn = @(btn,event) updateRTVFmodelParams(btn, obj);

    obj.theRTVFmodelTable.Layout.Row = 5;
    obj.theRTVFmodelTable.Layout.Column = 2;
    obj.theRTVFmodelTable.RowName = {};
    obj.theRTVFmodelTable.ColumnName = {'employed H1 horizontal cell index', 'grid half samples num, N (total samples: (2*N + 1)^2)'};
    obj.theRTVFmodelTable.Data = {};
    obj.theRTVFmodelTable.FontSize = 18;
    obj.theRTVFmodelTable.BackgroundColor = [0.4 0.4 0.4];
    addStyle(obj.theRTVFmodelTable, uistyle('HorizontalAlignment','center', 'FontColor', [0.3 1 0.8]))


    % The Action label
    theActionLabel.Layout.Row = 1;
    theActionLabel.Layout.Column = 1;
    theActionLabel.Text = "pipeline";
    theActionLabel.HorizontalAlignment = "center";
    theActionLabel.FontSize = 20;
    theActionLabel.FontColor = [0.1 0.8 0.9];
    theActionLabel.FontWeight = 'Bold';

    % The Action dropdown
    theActionDropdown.Layout.Row = 1;
    theActionDropdown.Layout.Column = 2;
    theActionDropdown.FontSize = 14;
    theActionDropdown.BackgroundColor = [0.3 0.8 0.95];
    theActionDropdown.Items = [ ...
        "compute: center-connected mRGC mosaic", ...
        "compute: all R2VFT objects", ...
        "inspect: single R2VFT object file", ...
        "inspect: all R2VFT objects file", ...
        "compute: refit R2VFT object(s)", ...
        "compute: manually replace a specific R2VFT object", ...
        "export: retinal cone pooling params for all fitted R2VFT objects", ...
        "compute: center-surround cone pooling kernels", ...
        "visualize: spatial RFs", ...
        "compute: frozen midgetRGCMosaic (with current center-surround cone pooling weights)", ...
        "validate: pre-compute LM-non-opponent cone mosaic STFs", ...
        "validate: compute LM-non-opponent RGC STFs (midgetRGCMosaic object handles the computation of cone mosaic responses)", ...
        "validate: compute LM-non-opponent RGC STFs from pre-computed LM-non-opponent cone mosaic STFs", ...
        "validate: fit LM-non-opponent STFs", ...
        "validate: visualize STF fits" ...
        ];

    % Current action
    obj.currentAction = theActionDropdown.Items{4};
    theActionDropdown.Value = obj.currentAction;
    theActionDropdown.ValueChangedFcn = @(src,event) dropDownActionChanged(src,event, obj);
    
    % The ExecuteActionButton
    theExecuteActionButton.Layout.Row = 2;
    theExecuteActionButton.Layout.Column = 2;
    theExecuteActionButton.Text = "commit selected pipeline";
    theExecuteActionButton.FontSize = 20;
    theExecuteActionButton.FontWeight = 'Bold';
    theExecuteActionButton.BackgroundColor = [0.4 0.4 0.4];
    theExecuteActionButton.FontColor = [0.1 0.8 0.9];
    theExecuteActionButton.ButtonPushedFcn = @(btn,event) executeButtonAction(btn, obj);


end

function updateRTVFmodelParams(btn, app)
end


function updateOpticsParams(btn, app)
end


function updateMosaicParams(btn, app)
end

function exitButtonAction(btn, app)
    app.mainView.delete();
end