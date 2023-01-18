classdef midgetRGCMosaicGeneratorApp < handle
    properties  (GetAccess=public, SetAccess=private)
        % GUI components
        mainView;

        % Current action to perform
        currentAction;

        % State (mosaic and optics params)
        simulation;
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
    % Load the various params structs
    [obj.simulation.mosaicCenterParams, ...
     obj.simulation.mosaicSurroundParams, ...
     obj.simulation.opticsParams] = midgetRGCMosaicGenerator.generateMosaicAndOpticsParamStructs();
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

            RTVobjIndicesToBeComputed = [9:22]
            pause
            
            midgetRGCMosaicGenerator.generateR2VFTobjects(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.mosaicSurroundParams, ...
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

            % Re-generate R2VFT objects at a specific position
            midgetRGCMosaicGenerator.generateR2VFTobjects(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.mosaicSurroundParams, ...
                app.simulation.opticsParams, ...
                'updateRTVFobjectAtPosition', targetPosition, ...
                'updateRTVFobjectWithCenterConesNum', targetRFcenterConesNum);

        case "compute: manually replace a specific R2VFT object"
            midgetRGCMosaicGenerator.replaceSpecificR2VFTobject()

        case "export: retinal cone pooling params for all fitted R2VFT objects"
            midgetRGCMosaicInspector.exportRetinalConePoolingParamsForAllFittedRTVFTobjects();

        case "compute: center-surround cone pooling kernels"
            midgetRGCMosaicGenerator.generateCenterSurroundConePoolingKernels(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.mosaicSurroundParams);

        case "visualize: spatial RFs"
            % Ask user how many RGCs to visualize
            maxRGCsNum = input('How many RGC RFs to visualize? ([Hit enter for all): ');

            midgetRGCMosaicInspector.visualizeSpatialRFmaps(...
                app.simulation.mosaicCenterParams, ...
                maxRGCsNum);
         
        case "compute: frozen mRGC mosaic"
            midgetRGCMosaicGenerator.freezeMosaic(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.mosaicSurroundParams);

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
                 app.simulation.mosaicSurroundParams, ...
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
                 app.simulation.mosaicSurroundParams, ...
                 useParfor);

        case "validate: compute LM-non-opponent RGC STFs from pre-computed LM-non-opponent cone mosaic STFs"
            useParfor = true;
            midgetRGCMosaicInspector.computeMosaicLMnonOpponentSTFsFromConeMosaicLMnonOpponentSTFs(...
                 app.simulation.mosaicCenterParams, ...
                 app.simulation.mosaicSurroundParams, ...
                 useParfor);

        case "validate: fit LM-non-opponent STFs"
            % Ask user how many RGCs to analyze
            maxRGCsNum = input('How many RGCs to fit? ([Hit enter for all): ');

            midgetRGCMosaicInspector.fitMosaicSTFs(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.mosaicSurroundParams, ...
                maxRGCsNum);


        otherwise
            error('Unknown action: ''%s''.', app.currentAction);
    end
end

function dropDownActionChanged(src,event, app)
    app.currentAction = event.Value;
end

function generateGUI(obj)

    % Create figure window
    obj.mainView = uifigure('Position', [30 500 1200 200], ...
        'WindowStyle', 'AlwaysOnTop', ...
        'Scrollable', 'on', ...
        'Resize', 'off');
    obj.mainView.Name = "Midget RGC Mosaic Generator & Inspector";

    % Manage app layout
    layoutRows = 1;
    layoutCols = 4;
    gl = uigridlayout(obj.mainView,[layoutRows layoutCols]);
    gl.RowHeight = {'1x', '4x'};
    gl.ColumnWidth = {'fit','6x', '1x', '1x'};

    theActionLabel  = uilabel(gl);
    theActionDropdown = uidropdown(gl);
    theExecuteActionButton = uibutton(gl);
    theExitButton = uibutton(gl);

    % The Action label
    theActionLabel.Layout.Row = 1;
    theActionLabel.Layout.Column = 1;
    theActionLabel.Text = "Actions:";
    theActionLabel.FontSize = 14;

    % The Action dropdown
    theActionDropdown.Layout.Row = 1;
    theActionDropdown.Layout.Column = 2;
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
        "compute: frozen mRGC mosaic", ...
        "validate: pre-compute LM-non-opponent cone mosaic STFs", ...
        "validate: compute LM-non-opponent RGC STFs (midgetRGCMosaic object handles the computation of cone mosaic responses)", ...
        "validate: compute LM-non-opponent RGC STFs from pre-computed LM-non-opponent cone mosaic STFs", ...
        "validate: fit LM-non-opponent STFs" ...
        ];

    theActionDropdown.FontSize = 14;
    % Current action
    obj.currentAction = theActionDropdown.Items{4};
    theActionDropdown.Value = obj.currentAction;
    theActionDropdown.ValueChangedFcn = @(src,event) dropDownActionChanged(src,event, obj);
    
    % The ExecuteActionButton
    theExecuteActionButton.Layout.Row = 1;
    theExecuteActionButton.Layout.Column = 3;
    theExecuteActionButton.Text = "GO !";
    theExecuteActionButton.FontSize = 20;
    theExecuteActionButton.BackgroundColor = [0.3 0.3 0.3];
    theExecuteActionButton.FontColor = [0.1 0.8 0.9];
    theExecuteActionButton.ButtonPushedFcn = @(btn,event) executeButtonAction(btn, obj);

    % The ExitButton
    theExitButton.Layout.Row = 1;
    theExitButton.Layout.Column = 4;
    theExitButton.Text = "Quit";
    theExitButton.FontSize = 20;
    theExitButton.BackgroundColor = [0.3 0.3 0.3];
    theExitButton.FontColor = [0.9 0.8 0.0];
    theExitButton.ButtonPushedFcn = @(btn,event) exitButtonAction(btn, obj);
end

function exitButtonAction(btn, app)
    app.mainView.delete();
end