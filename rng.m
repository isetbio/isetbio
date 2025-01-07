function settings = rng(arg1,arg2)

    global rngTrackingInfo

    dialogHeader = '';
    theUIFigure = [];

    if (isfield(rngTrackingInfo, 'scenarioBeingRun'))

        % Update rng call no
        rngTrackingInfo.callNo = rngTrackingInfo.callNo + 1;

        scenarioBeingRun = rngTrackingInfo.scenarioBeingRun;
        codePathBeingTested = rngTrackingInfo.rngCodePathToRun;
        rngCallNo = rngTrackingInfo.callNo;
        

        switch (nargin)
            case 0
                dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with no arguments. \n', scenarioBeingRun, codePathBeingTested, rngCallNo);
    
            case 1
                if (ischar(arg1))
                    dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with a single argument: ''%s''.', scenarioBeingRun, codePathBeingTested, rngCallNo, arg1);
    
                elseif (isnumeric(arg1))
                    dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with a single argument: %g', scenarioBeingRun, codePathBeingTested, rngCallNo, arg1);
    
                elseif (isstruct(arg1))
                    dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with a struct with the following fields:', scenarioBeingRun, codePathBeingTested, rngCallNo);
                    fNames = fieldnames(arg1);
                    for i = 1:numel(fNames)
                        if (ischar(arg1.(fNames{i})))
                            dialogHeader = sprintf('%s\n%s: ''%s''', dialogHeader, fNames{i}, arg1.(fNames{i}));
                        elseif (numel(arg1.(fNames{i})) == 1)
                            dialogHeader = sprintf('%s\n%s: %g', dialogHeader, fNames{i}, arg1.(fNames{i}));
                        else
                            dialogHeader = sprintf('%s\n%s, containing %d elements', dialogHeader, fNames{i}, numel(arg1.(fNames{i})));
                        end
                    end
    
                else
                    dialogHeader = sprintf('Running scenario: ''%s'' with path: ''%s''.\n\nrng() call (#%d) with an argument that is neither a char, a numeric, or a struct', scenarioBeingRun, codePathBeingTested);
                    arg1
                    class(arg1)
                end
    
            case 2
    
                if (isnumeric(arg1)&&isnumeric(arg2))
                    dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with 2 arguments: %g %g', scenarioBeingRun, codePathBeigTested, rngCallNo, arg1, arg2);
                else
                    dialogHeader = sprintf('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with 2 arguments that are not numeric', scenarioBeingRun, codePathBeingTested, rngCallNo);
                    arg1
                    arg2
                    disp('classes of arg1, arg2')
                    class(arg1)
                    class(arg2)
                end
    
            otherwise
                error('Running scenario: ''%s''\nCode path tested: ''%s''.\n\nrng() call (#%d) with %d arguments. Dont know how to handle this.', scenarioBeingRun, codePathBeingTested, rngCallNo, nargin)
    
        end % switch (nargin)


        % Display the calling stack in a separate window
        theUIFigure = rngTrackingInfo.callingStackUIFigure;

        if (~isempty(theUIFigure))
            % Assemble calling stack info
            theRootDir = strrep(isetRootPath, 'isetcam', '');
            stack = dbstack('-completenames');
            
            callingStackInfo = {};
            theCallingFunctionNamesTableColumn = {};
            theFilenamesTableColumn = {};
            theLineNumbersTableColumn = [];
            theCallingOrder = {};
        
            for i = numel(stack):-1:1
                s = stack(i);
                theFileName = strrep(strrep(s.file, theRootDir,''), 'isetbio', '');
                theLineNo = s.line;
                theCallingFunctionName = s.name;
                callingStackInfo{numel(callingStackInfo)+1} = sprintf('%d) %s->%s (%d)', ...
                    numel(stack)-i+1, theFileName, theCallingFunctionName, theLineNo);
                theFilenamesTableColumn{numel(stack)-i+1} = theFileName;
                theCallingFunctionNamesTableColumn{numel(stack)-i+1} = theCallingFunctionName;
                theLineNumbersTableColumn{numel(stack)-i+1} = theLineNo;
                theCallingOrder{numel(stack)-i+1} = sprintf('%d)', numel(stack)-i+1);
            end
    
   
            T = table(...
                string(theCallingOrder'), ...
                string(theFilenamesTableColumn'), ...
                string(theCallingFunctionNamesTableColumn'),...
                 theLineNumbersTableColumn');
        
            T.Properties.VariableNames = ["Call no", "Filename", "Function name", "Line no"];
             
            uit = uitable(theUIFigure, ...
                'Data',T, ...
                'Position',[20 20 1800 300], ...
                'FontSize', 16);
            set (uit,'ColumnWidth', {80,700,600,100});
        end % if (~isempty(theUIFigure))

    end


    % Display the current call to the rng
    % Get the default UI font size
    defaultUIcontrolFontSize = get(groot,'defaultUicontrolFontSize');
    defaultUIcontrolFontName = get(groot, 'defaultUIcontrolFontName');
    
    % Set bigger font
    set(groot, 'defaultUicontrolFontSize', 16);
    set(groot, 'defaultUIcontrolFontName', 'Menlo');

    % Present dialog
    options.Default = ' proceed with execution ';
    options.Interpreter = 'none';

    % This will hang here until we enter a choice
    response = questdlg(dialogHeader, ...
	    'INTERCEPTING RNG() CALL PATH', ...
	    ' proceed with execution ',  ' interrupt execution here ', ...
        options);

    % Restore previous font size
    set(groot,'defaultUicontrolFontSize',defaultUIcontrolFontSize);
    set(groot,'defaultUicontrolFontName',defaultUIcontrolFontName);

    % Save current directory
    oldDir = pwd;



    % Handle user's response
    switch response
        case ' proceed with execution '

            % Move to built-in rng directory so we can execute it
            cd('/Applications/MATLAB_R2024b.app/toolbox/matlab/randfun');

            % Dispatch rng call
            switch (nargin)
                case 0
                    settings = feval('rng');
                case 1
                    settings = feval('rng',arg1);
                case 2
                    settings = feval('rng',arg1,arg2);
            end

            % Restore original directory
            cd(oldDir);

        case ' interrupt execution here '

            % Restore original directory
            cd(oldDir);

            error('Interrupted execution at current point');
    end % switch response

end

