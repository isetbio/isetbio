classdef midgetRGCMosaicGeneratorApp < handle


    properties  (GetAccess=public, SetAccess=private)
        % GUI components
        mainView;

        % Current action to perform
        currentPipeline;

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


function importParams(btn, app)

    suggestedParamsFileDirectory = fullfile(isetRootPath, 'ganglioncells');

    [fileName, filePath] = uigetfile({'*.mat'} ,...
                        'Select a params file', suggestedParamsFileDirectory);

    theUserSelectedParamsFile = fullfile(filePath, fileName);
    load(theUserSelectedParamsFile, ...
             'mosaicCenterParams', ...
             'rfModelParams', ...
             'opticsParams');

    app.simulation.mosaicCenterParams = mosaicCenterParams;
    app.simulation.rfModelParams = rfModelParams;
    app.simulation.opticsParams = opticsParams;

    updateParamsTables(app);
end

function exportParams(btn, app)
    % Retrieve current params
    mosaicCenterParams = app.simulation.mosaicCenterParams;
    rfModelParams = app.simulation.rfModelParams;
    opticsParams = app.simulation.opticsParams;

    suggestedParamsFileName = 'params.mat';
    suggestedParamsFileDirectory = fullfile(isetRootPath, 'ganglioncells');

    [fileName, filePath] = uiputfile(...
        suggestedParamsFileName , 'Specify params filename', ...
        suggestedParamsFileDirectory);
    if ((isequal(fileName,0)) || (isequal(filePath,0)))
        disp('Not exporting anything\n');
        return;
    end

    theUserSelectedParamsFile = fullfile(filePath, fileName);
    save(theUserSelectedParamsFile, ...
             'mosaicCenterParams', ...
             'rfModelParams', ...
             'opticsParams');
end


function initializeState(obj)

     theParamsFile = midgetRGCMosaicGenerator.paramsFile();

     if (isfile(theParamsFile))
         fprintf('Loading previously used params...\n');
         load(theParamsFile, ...
             'mosaicCenterParams', ...
             'rfModelParams', ...
             'opticsParams');
         obj.simulation.mosaicCenterParams = mosaicCenterParams;
         obj.simulation.rfModelParams  = rfModelParams;
         obj.simulation.opticsParams = opticsParams;
     else
        % Generate the default params structs
        fprintf('Loading default params...\n');
        [obj.simulation.mosaicCenterParams, ...
         obj.simulation.rfModelParams, ...
         obj.simulation.opticsParams] = midgetRGCMosaicGenerator.generateDefaultMosaicAndOpticsParamStructs();
     end

     updateParamsTables(obj);
end


function updateParamsTables(obj)
     % Update theMosaicGeometryTable
     d = {};
     d{1} = obj.simulation.mosaicCenterParams.positionDegs(1);
     d{2} = obj.simulation.mosaicCenterParams.positionDegs(2);
     d{3} = obj.simulation.mosaicCenterParams.sizeDegs(1);
     d{4} = obj.simulation.mosaicCenterParams.sizeDegs(2);
     set(obj.theMosaicGeometryTable,'Data',d)

     % Update theOpticsTable
     d = {};
     d{1} = obj.simulation.opticsParams.ZernikeDataBase;
     d{2} = obj.simulation.opticsParams.subjectRankOrder;
     d{3} = obj.simulation.opticsParams.pupilDiameterMM;
     set(obj.theOpticsTable, 'Data', d);

     % Update rfModelParamsTable
     d = {};
     d{1} = obj.simulation.rfModelParams.H1cellIndex;
     d{2} = obj.simulation.rfModelParams.eccentricitySamplingGridHalfSamplesNum;
     set(obj.theRTVFmodelTable, 'Data', d)
     drawnow
end

function executePipelineAction(btn, app)

    switch app.currentPipeline
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
            error('Unknown action: ''%s''.', app.currentPipeline);
    end
end

function dropDownPipelineChanged(src,event, app)
    app.currentPipeline = event.Value;
end

function generateGUI(obj)

    % Create figure window
    obj.mainView = uifigure('Position', [30 1500 1000 420], ...
        'WindowStyle', 'AlwaysOnTop', ...
        'Scrollable', 'on', ...
        'Color', [0.3 0.3 0.3], ...
        'Resize', 'off');

    obj.mainView.Name = "midgetRGCMosaic generator & inspector";

    % Manage app layout
    layoutRows = 6;
    layoutCols = 2;
    gl = uigridlayout(obj.mainView,[layoutRows layoutCols]);
    gl.BackgroundColor = [0.3 0.3 0.3];
    gl.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x'};
    gl.ColumnWidth = {'fit','6x'};

    thePipelineLabel  = uilabel(gl);
    theMosaicGeometryButton  = uibutton(gl);
    obj.theMosaicGeometryTable = uitable(gl);
    obj.theOpticsTable = uitable(gl);
    obj.theRTVFmodelTable = uitable(gl);
    theOpticsButton = uibutton(gl);
    theRTVFmodelButton = uibutton(gl);
    thePipelineDropdown = uidropdown(gl);
    theExecutePipelineButton = uibutton(gl);
    theExportParamsButton = uibutton(gl);
    theImportParamsButton = uibutton(gl);

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
    obj.theMosaicGeometryTable.ColumnName = {'POSITION,X (DEGS)', 'POSITION, Y (degs)', 'SIZE, X (degs)', 'SIZE, Y(degs)'};
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
    obj.theOpticsTable.ColumnName = {'ZERNIKE DATABASE', 'SUBJECT RANK', 'PUPIL DIAMETER (MM)'};
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
    obj.theRTVFmodelTable.ColumnName = {'EMPLOYED H1 HORIZONTAL CELL INDEX (SURROUND)', 'SPATIAL GRID HALF SAMPLES NUM, N (TOTAL: (2*N + 1)^2)'};
    obj.theRTVFmodelTable.Data = {};
    obj.theRTVFmodelTable.FontSize = 18;
    obj.theRTVFmodelTable.BackgroundColor = [0.4 0.4 0.4];
    addStyle(obj.theRTVFmodelTable, uistyle('HorizontalAlignment','center', 'FontColor', [0.3 1 0.8]))


    % The Pipelinelabel
    thePipelineLabel.Layout.Row = 1;
    thePipelineLabel.Layout.Column = 1;
    thePipelineLabel.Text = "pipeline";
    thePipelineLabel.HorizontalAlignment = "center";
    thePipelineLabel.FontSize = 20;
    thePipelineLabel.FontColor = [1 0.8 0.2];
    thePipelineLabel.FontWeight = 'Bold';

    % The Pipeline dropdown
    thePipelineDropdown.Layout.Row = 1;
    thePipelineDropdown.Layout.Column = 2;
    thePipelineDropdown.FontSize = 14;
    thePipelineDropdown.BackgroundColor = [1 0.9 0.6];
    thePipelineDropdown.Items = [ ...
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
    obj.currentPipeline = thePipelineDropdown.Items{4};
    thePipelineDropdown.Value = obj.currentPipeline;
    thePipelineDropdown.ValueChangedFcn = @(src,event) dropDownPipelineChanged(src,event, obj);
    
    % The ExecutePipelineButton
    theExecutePipelineButton.Layout.Row = 2;
    theExecutePipelineButton.Layout.Column = 2;
    theExecutePipelineButton.Text = "commit selected pipeline";
    theExecutePipelineButton.FontSize = 20;
    theExecutePipelineButton.FontWeight = 'Bold';
    theExecutePipelineButton.BackgroundColor = [0.4 0.4 0.4];
    theExecutePipelineButton.FontColor = [1 0.8 0.2];
    theExecutePipelineButton.ButtonPushedFcn = @(btn,event) executePipelineAction(btn, obj);


    % TheExportParamsButton
    theExportParamsButton.Layout.Row = 6;
    theExportParamsButton.Layout.Column = 1;
    theExportParamsButton.Text = "export params";
    theExportParamsButton.FontSize = 20;
    theExportParamsButton.FontWeight = 'Bold';
    theExportParamsButton.BackgroundColor = [0.4 0.4 0.4];
    theExportParamsButton.FontColor = [0.3 1 0.8];
    theExportParamsButton.ButtonPushedFcn = @(btn,event) exportParams(btn, obj);

    % TheImportParamsButton
    theImportParamsButton.Layout.Row = 6;
    theImportParamsButton.Layout.Column = 2;
    theImportParamsButton.Text = "import params";
    theImportParamsButton.FontSize = 20;
    theImportParamsButton.FontWeight = 'Bold';
    theImportParamsButton.BackgroundColor = [0.4 0.4 0.4];
    theImportParamsButton.FontColor = [0.3 1 0.8];
    theImportParamsButton.ButtonPushedFcn = @(btn,event) importParams(btn, obj);


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