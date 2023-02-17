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

    midgetRGCMosaicInspector.say('Choose the params file to import');
    [fileName, filePath] = uigetfile({'*.mat'} ,...
                        'Choose the params file to import', suggestedParamsFileDirectory);

    if ((isequal(fileName,0)) || (isequal(filePath,0)))
        disp('Import cancelled');
        return;
    end

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

    midgetRGCMosaicInspector.say('Choose exported params filename');
    
    [fileName, filePath] = uiputfile(...
        suggestedParamsFileName , 'Choose exported params filename', ...
        suggestedParamsFileDirectory);
    if ((isequal(fileName,0)) || (isequal(filePath,0)))
        disp('Export cancelled');
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
        case "GENERATION: center-connected mRGC mosaic"
            midgetRGCMosaicGenerator.generateCenterConnectedMosaic(...
                app.simulation.mosaicCenterParams);


        case "OPTIMIZATION: fit all locations R2VFT objects"

            RTVobjIndicesToBeComputed = 'all';
            midgetRGCMosaicGenerator.fitMultiFocalRTVF(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                'RTVobjIndicesToBeComputed', RTVobjIndicesToBeComputed);


        case "VISUALIZATION: fits in single location R2VFT object file"
            midgetRGCMosaicInspector.quicklyInspectSingleRTVFobjectFile();

        case "VISUALIZATION: fits in all locations R2VFT objects file"
            midgetRGCMosaicInspector.quicklyInspectAllRTVFobjectsFile();

        case "OPTIMIZATION: refit single location R2VFT object"
            % Ask user which R2VFT objects to refit 
            % (single position,  multiple center cones num)

            singleRTVobjIndexToBeComputed = [];
            while (isempty(singleRTVobjIndexToBeComputed))
                singleRTVobjIndexToBeComputed = input('\nEnter index of RTVF object to be refit: ');
            end

            % Ask user whether to compute the L-cone compute struct
            computeLconeCenterComputeStruct = queryUserForYesNoResponse('\nCompute the L-cone compute struct?');

            % Ask user whether to compute the M-cone compute struct
            computeMconeCenterComputeStruct = queryUserForYesNoResponse('\nCompute the M-cone compute struct?');

            multiStartsNumRetinalPooling  = [];
            while (isempty(multiStartsNumRetinalPooling))
                multiStartsNumRetinalPooling = input('\nNumber of multi-starts, e.g., 1, 2, ... : ');
            end

            midgetRGCMosaicGenerator.fitMultiFocalRTVF(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                'RTVobjIndicesToBeComputed', singleRTVobjIndexToBeComputed, ...
                'computeLconeCenterComputeStruct', computeLconeCenterComputeStruct, ...
                'computeMconeCenterComputeStruct', computeMconeCenterComputeStruct, ...
                'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling);

        case "EXPORT: optimized retinal cone pooling params for all fitted R2VFT objects"
            midgetRGCMosaicInspector.exportRetinalConePoolingParamsForAllFittedRTVFTobjects();

        case "GENERATION: center-surround cone pooling kernels"
            midgetRGCMosaicGenerator.generateCenterSurroundConePoolingKernels(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams)

        case "VISUALIZATION: multi-focal RTVF L- and M-cone weighted PSFs"
            midgetRGCMosaicInspector.visualizeFittedMultiFocalRTVF(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams);
            
        case "VISUALIZATION: spatial RFs along arbitrary meridians"
            % Ask user how many RGCs to visualize
            maxRGCsNum = input('How many RGC RFs to visualize? ([Hit enter for all): ');

            midgetRGCMosaicInspector.visualizeSpatialRFmaps(...
                app.simulation.mosaicCenterParams, ...
                maxRGCsNum);
         
        case "GENERATION: frozen midgetRGCMosaic (with current center-surround cone pooling weights)"
            midgetRGCMosaicGenerator.freezeMosaic(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams)

        case "VALIDATION: pre-compute cone mosaic STFs"
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

        case "VALIDATION: compute RGC STFs (midgetRGCMosaic object handles the computation of cone mosaic responses)"
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

        case "VALIDATION: compute RGC STFs from pre-computed cone mosaic STFs"
            useParfor = true;
            midgetRGCMosaicInspector.computeMosaicLMnonOpponentSTFsFromConeMosaicLMnonOpponentSTFs(...
                 app.simulation.mosaicCenterParams, ...
                 app.simulation.rfModelParams, ...
                 app.simulation.opticsParams, ...
                 useParfor);

        case "VALIDATION: fit STFs"
            % Ask user how many RGCs to analyze
            maxRGCsNum = input('How many RGCs to fit? ([Hit enter for all): ');

            midgetRGCMosaicInspector.fitMosaicSTFs(...
                app.simulation.mosaicCenterParams, ...
                app.simulation.rfModelParams, ...
                app.simulation.opticsParams, ...
                maxRGCsNum);


        case "VALIDATION: visualize STF fits across eccentricities"
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
    obj.mainView = uifigure('Position', [30 1500 1000 450], ...
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
    obj.theMosaicGeometryTable.ColumnName = {'POSITION,X (DEGS)', 'POSITION, Y (DEGS', 'SIZE, X (DEGS)', 'SIZE, Y(DEGS)'};
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
    theRTVFmodelButton.Text = "RTVF params";
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

    % The Pipeline dropdown
    thePipelineDropdown.Layout.Row = 1;
    thePipelineDropdown.Layout.Column = 2;
    thePipelineDropdown.FontSize = 14;
    thePipelineDropdown.BackgroundColor = [1 0.9 0.6];
    thePipelineDropdown.Items = [ ...
        "GENERATION: center-connected mRGC mosaic", ...
        "OPTIMIZATION: fit all locations R2VFT objects", ...
        "VISUALIZATION: fits in single location R2VFT object file", ...
        "VISUALIZATION: fits in all locations R2VFT objects file", ...
        "OPTIMIZATION: refit single location R2VFT object", ...
        "EXPORT: optimized retinal cone pooling params for all fitted R2VFT objects", ...
        "GENERATION: center-surround cone pooling kernels", ...
        "VISUALIZATION: multi-focal RTVF L- and M-cone weighted PSFs", ...
        "VISUALIZATION: spatial RFs along arbitrary meridians", ...
        "GENERATION: frozen midgetRGCMosaic (with current center-surround cone pooling weights)", ...
        "VALIDATION: pre-compute cone mosaic STFs", ...
        "VALIDATION: compute RGC STFs (midgetRGCMosaic object handles the computation of cone mosaic responses)", ...
        "VALIDATION: compute RGC STFs from pre-computed cone mosaic STFs", ...
        "VALIDATION: fit STFs", ...
        "VALIDATION: visualize STF fits across eccentricities" ...
        ];
    % Current action
    obj.currentPipeline = thePipelineDropdown.Items{4};
    thePipelineDropdown.Value = obj.currentPipeline;
    thePipelineDropdown.ValueChangedFcn = @(src,event) dropDownPipelineChanged(src,event, obj);
    
    
    % The ExecutePipelineButton
    theExecutePipelineButton.Layout.Row = 1;
    theExecutePipelineButton.Layout.Column = 1;
    theExecutePipelineButton.Text = "commit";
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
    theImportParamsButton.Layout.Row = 2;
    theImportParamsButton.Layout.Column = 1;
    theImportParamsButton.Text = "import params";
    theImportParamsButton.FontSize = 20;
    theImportParamsButton.FontWeight = 'Bold';
    theImportParamsButton.BackgroundColor = [0.4 0.4 0.4];
    theImportParamsButton.FontColor = [0.3 1 0.8];
    theImportParamsButton.ButtonPushedFcn = @(btn,event) importParams(btn, obj);


end

function updateRTVFmodelParams(btn, app)
    
   app.simulation.rfModelParams.H1cellIndex = queryUserForParamValue(...
        'employed H1 cell index', [1 2 3 4], app.simulation.rfModelParams.H1cellIndex);
    updateParamsTables(app);

   app.simulation.rfModelParams.eccentricitySamplingGridHalfSamplesNum = queryUserForParamValue(...
        'spatial grid half-samples num', [], app.simulation.rfModelParams.eccentricitySamplingGridHalfSamplesNum );
   updateParamsTables(app);
   

end


function updateOpticsParams(btn, app)
    app.simulation.opticsParams.ZernikeDataBase = queryUserForParamValue(...
        'Zernike database', {'Polans2015', 'Artal2012'}, app.simulation.opticsParams.ZernikeDataBase);
    updateParamsTables(app);

    app.simulation.opticsParams.subjectRankOrder = queryUserForParamValue(...
        'subject rank', [], app.simulation.opticsParams.subjectRankOrder);
    updateParamsTables(app);

    app.simulation.opticsParams.pupilDiameterMM = queryUserForParamValue(...
        'pupil diameter (mm)', [], app.simulation.opticsParams.pupilDiameterMM);
    updateParamsTables(app);
end


function updateMosaicParams(btn, app)

    app.simulation.mosaicCenterParams.positionDegs(1) = queryUserForParamValue(...
        'position, x (degs)', [], app.simulation.mosaicCenterParams.positionDegs(1));
    updateParamsTables(app);

    app.simulation.mosaicCenterParams.positionDegs(2) = queryUserForParamValue(...
        'position, y (degs)', [], app.simulation.mosaicCenterParams.positionDegs(2));
    updateParamsTables(app);

    app.simulation.mosaicCenterParams.sizeDegs(1) = queryUserForParamValue(...
        'size, x (degs)', [], app.simulation.mosaicCenterParams.sizeDegs(1));
    updateParamsTables(app);

    app.simulation.mosaicCenterParams.sizeDegs(2) = queryUserForParamValue(...
        'size, y (degs)', [], app.simulation.mosaicCenterParams.sizeDegs(2));
    updateParamsTables(app);
end

function val = queryUserForParamValue(paramName, validParamValues, oldVal)
    val = oldVal;
    if (~isempty(validParamValues))
        fprintf('Valid values for ''%s'': \n', paramName);
        validParamValues
    end
    if (ischar(oldVal))
        newVal = input(sprintf('Update ''%s'' (current: ''%s'')? Enter to keep current value:', paramName, oldVal), 's');
    else
        newVal = input(sprintf('Update ''%s'' (current: %g) ? Enter to keep current value: ', paramName, oldVal));
    end

    if (~isempty(newVal))
        if ((~isempty(validParamValues)) && (ismember(newVal, validParamValues))) || (isempty(validParamValues))
            val = newVal;
        end
    end
end

function val = queryUserForYesNoResponse(message)
    notValidResponse = true;
    while (notValidResponse)
        message = sprintf('%s [y=YES, n=NO]: ', message);
        txt = lower(input(message, 's'));
        if (strcmp(txt, 'y')) || (strcmp(txt, 'n'))
            notValidResponse = false;
        end
    end   
    if (strcmp(txt, 'y'))
        val = true;
        fprintf('Will do.\n');
    else
        val = false;
        fprintf('Will SKIP\n');
   end
end


function exitButtonAction(btn, app)
    app.mainView.delete();
end