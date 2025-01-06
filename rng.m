function settings = rng(arg1,arg2)

    global scenarioBeingRun


    switch (nargin)
        case 0
            dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with no arguments. \n', scenarioBeingRun);

        case 1
            if (ischar(arg1))
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with a single argument: ''%s''.', scenarioBeingRun, arg1);

            elseif (isnumeric(arg1))
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with a single argument: %g', scenarioBeingRun, arg1);

            elseif (isstruct(arg1))
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with a struct with the following fields:', scenarioBeingRun);
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
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with an argument that is neither a char, a numeric, or a struct', scenarioBeingRun);
                arg1
                class(arg1)
            end

        case 2

            if (isnumeric(arg1)&&isnumeric(arg2))
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with 2 arguments: %g %g', scenarioBeingRun, arg1, arg2);
            else
                dialogHeader = sprintf('Running scenario: ''%s''.\n>> rng called with 2 arguments that are not numeric', scenarioBeingRun);
                arg1
                arg2
                disp('classes of arg1, arg2')
                class(arg1)
                class(arg2)
            end

        otherwise
            error('Running scenario: ''%s''.\n>> rng called with %d arguments. Dont know how to handle this.', scenarioBeingRun, nargin)

    end % switch (nargin)


    % Get the default UI font size
    defaultUIcontrolFontSize = get(groot,'defaultUicontrolFontSize');

    % Set big font
    set(groot,'defaultUicontrolFontSize', 18);


    % Present dialog
    options.Default = ' proceed with execution ';
    options.Interpreter = 'none';
    response = questdlg(dialogHeader, ...
	    'INTEREPTING RNG() CALL', ...
	    ' proceed with execution ',  ' interrupt execution to see the code path ', ...
        options);

    % Restore previous font size
    set(groot,'defaultUicontrolFontSize',defaultUIcontrolFontSize);

    % Handle user's response
    switch response
        case ' proceed with execution '

            % Save current directory
            oldDir = pwd;

            % Move to built in rng directory so we can execute it
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

        case ' interrupt execution to see the code path '
            error('Interrupted execution to see the code path');
    end % switch response

end


