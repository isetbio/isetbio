function outContainer = combineContainers(theCellArrayOfContainers)
% Consolidate cell array of containers into a container of a cell array.
% 
% Syntax:
%    outContainer = combineContainers(theCellArrayOfContainers)
%
% Description:
%    Function that takes a cell array of containers and returns a single
%    container whose contents are a cell array of what was in each of
%    the passed containers.  This operates over all the flags of the
%    passed containers, and check that the flags are matched across the
%    containers passed.
%
% Inputs:
%    theCellArrayOfContainers      - The passed cell array of containers
%
% Outputs:
%    outContainer                  - The consolidated container of cell
%                                    arrays.
%
% Optional key/value pairs:
%    None
%
% See also: combineContainersMat
%

% History:
%    12/12/21  ncp  Wrote it.
%    12/13/21  dhb  More comments.

% Examples:
%{
    % Keys common to all containers
    monthNames = {'Jan','Feb','Mar','Apr', 'May', 'June'};
    
    % Values for container monthDays
    days = [31 28 31 30 31 30];
    monthDays = containers.Map(monthNames,days);
    
    % Values for container monthIndex
    indices = 1:numel(monthNames);
    monthIndex = containers.Map(monthNames, indices);
    
    % Values for container monthGreekNames
    greekNames = {'Ianouarios', 'Februarios', 'Martios', 'Aprilios', 'Maios', 'Iounios'};
    monthGreekNames = containers.Map(monthNames, greekNames);
    
    % Construct cell array with all containers
    theCellArrayOfContainers{1} = monthDays;
    theCellArrayOfContainers{2} = monthIndex;
    theCellArrayOfContainers{3} = monthGreekNames;
    
    % Debugging
    makeItFail = [false false];
    if (makeItFail(1))
        % Make it fail (#1)
        theCellArrayOfContainers{4} = containers.Map(monthNames{1:3}, ['a', 'b', 'c']);
    end
    
    if (makeItFail(2))
        % Make it fail (#2)
        theCellArrayOfContainers{4} = containers.Map({'Jan','Feb','Mar','Apr', 'May', 'June', 'July'}, 1:7);
    end
    
    % Combine the containers
    outContainer = combineContainers(theCellArrayOfContainers);
    
    % Test that it worked. Note that the keys are not ordered in the initial order. 
    % That's a feature of containers.
    theKeys = keys(outContainer);
    for iKey = 1:numel(theKeys)
        % display the contents of the outContainer for each key
        theKey = theKeys{iKey};
        fprintf('Content of all containers for key: ''%s''', theKey);
        outContainer(theKey)
    end
%}

    % Ensure that all containers have the same flags
    validateContainerKeys(theCellArrayOfContainers);
    
    % Initialize output container
    outContainer = containers.Map();
    
    % Fill the entries of the outContainer
    theKeys = keys(theCellArrayOfContainers{1});
    for iEntry = 1:numel(theKeys)
        theKey = theKeys{iEntry};
        theKeyValuesAcrossAllContainers = cell(1,numel(theCellArrayOfContainers));
        for iContainer = 1:numel(theCellArrayOfContainers)
            theKeyValuesAcrossAllContainers{iContainer} = theCellArrayOfContainers{iContainer}(theKey);
        end
        outContainer(theKey) = theKeyValuesAcrossAllContainers;
    end
end

function validateContainerKeys(theCellArrayOfContainers)
% Ensure that all containers have the same flags

    for iContainer = 1:numel(theCellArrayOfContainers)
        theContainerKeys = keys(theCellArrayOfContainers{iContainer});
        if (iContainer == 1)
            theCommonKeys = theContainerKeys;
        else
            % Check that all keys in this container are also in theCommonKeys
            for iKey = 1:numel(theContainerKeys)
                assert(ismember(theContainerKeys{iKey}, theCommonKeys), 'key ''%s'' in the %d container does not exist in the first container', ...
                    theContainerKeys{iKey}, iContainer);
            end
            % Check that all keys in theCommonKeys are also keys in this container
            for iKey = 1:numel(theCommonKeys)
                assert(ismember(theCommonKeys{iKey}, theContainerKeys), 'key ''%s'' in the first container does not exist in the %d container', ...
                    theCommonKeys{iKey}, iContainer);
            end
        end
    end
end
