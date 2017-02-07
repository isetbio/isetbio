function result = RecursivelyCompareStructs(struct1Name, struct1, struct2Name, struct2, varargin)
% result = RecursivelyCompareStructsTests(struct1Name, struct1, struct2Name, struct2, defaultTolerance, customTolerances)
%
% Method to compare nested structs with arbitrary internal organization
% using a defaultTolerance as well as an optional customTolerance for any desired field.
%
% Simple usage:
% s1 = struct('a', 1);
% s2 = struct('a', 2);
% result = RecursivelyCompareStructs('s1', s1, 's2', s2)

% For more extreme usage examples see RecursivelyCompareStructsTests.
%
% 2/7/17  NPC  Wrote it.

%% Parse input
p = inputParser;
p.addParameter('defaultTolerance', 0, @isnumeric);
p.addParameter('customTolerances', containers.Map());
p.parse(varargin{:});
defaultTolerance = p.Results.defaultTolerance;
customTolerances = p.Results.customTolerances;

% Go
[result, resultOperations, ~, ~, ~] = doRecursion(...
    struct1Name, struct1, struct2Name, struct2, ...
    defaultTolerance, customTolerances, {}, {}, {}, ...
    [], [], struct2Name);

% Print the results
if (isempty(result))
    %fprintf('\n<strong>Structs comparison results: structs are identical.</strong>\n');
else
    fprintf('\n<strong>Structs comparison results: structs are NOT identical. </strong>\n');
    for k = 1:numel(result)
        fprintf('\t[%d]. %s\n', k, result{k});
    end
end

% More printing
if (1==2)
    for k = 1:numel(toleranceFieldPaths)
        fprintf('toleranceFieldPath{%d}: ''%s''\n', k, toleranceFieldPaths{k});
    end
    fprintf('\n');
    for k = 1:numel(flatStruct1FieldNames)
        fprintf('flatStruct1FieldNames{%d}: ''%s''\n', k, flatStruct1FieldNames{k});
    end

    fprintf('\n');
    for k = 1:numel(resultOperations)
        fprintf('resultOperation{%d}: ''%s''\n', k, resultOperations{k});
    end
end

end

function [result, resultOperations, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames] = doRecursion(...
    struct1Name, struct1, struct2Name, struct2, ...
    tolerance, customTolerances, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames, ...
     oldResult, oldResultOperations, topLevelStruct2Name)

    result = oldResult;
    resultOperations = oldResultOperations;

    % Check that structs have mathcing lengths
    assertStructsHaveMatchingLength(struct1Name, struct1, struct2Name, struct2);

    % OK, inputs are good structs so lets continue with their fields
    theStruct1Name = struct1Name;
    theStruct2Name = struct2Name;

    for structIndex = 1:numel(struct1)

        if (numel(struct1) > 1)
            struct1FullName = sprintf('%s(%d)', theStruct1Name, structIndex);
            struct2FullName = sprintf('%s(%d)', theStruct2Name, structIndex);
        else
            struct1FullName = sprintf('%s', theStruct1Name);
            struct2FullName = sprintf('%s', theStruct2Name);
        end
        struct1Name = sprintf('%s', theStruct1Name);
        struct2Name = sprintf('%s', theStruct2Name);
        if (strcmp(struct2Name, topLevelStruct2Name))
            structToleranceName = strrep(struct2Name, topLevelStruct2Name, '');
        else
            structToleranceName = struct2Name;
        end
        struct1FieldNames = sort(fieldnames(struct1));
        struct2FieldNames = sort(fieldnames(struct2));

        for fieldIndex = 1:numel(struct1FieldNames)

            fullFieldName1 = sprintf('%s.%s', struct1FullName,struct1FieldNames{fieldIndex});
            fullFieldName2 = sprintf('%s.%s', struct2FullName,struct2FieldNames{fieldIndex});

            field1 = [];
            field2 = [];
            
            eval(sprintf('field1 = struct1(structIndex).%s;', struct1FieldNames{fieldIndex}));
            eval(sprintf('field2 = struct2(structIndex).%s;', struct2FieldNames{fieldIndex}));

            % Update list of flat struct field names
            flatStruct1FieldNames{numel(flatStruct1FieldNames)+1} = fullFieldName1;
            flatStruct2FieldNames{numel(flatStruct2FieldNames)+1} = fullFieldName2;

            % Update toleranceFieldPaths
            if (structIndex == 1)
                idx = strfind(structToleranceName, '.');
                if (~isempty(idx))
                    structToleranceName = structToleranceName(idx(1)+1:end);
                end
                if (isempty(structToleranceName))
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s', struct1FieldNames{fieldIndex});
                else
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s.%s', structToleranceName,struct1FieldNames{fieldIndex});
                end
            end

            % Fields are structs
            if (isstruct(field1)) && (isstruct(field2))
                % Go into recursion
                [result, resultOperations, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames] = doRecursion(...
                    fullFieldName1, field1, fullFieldName2, field2, ...
                    tolerance, customTolerances, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames, ...
                    result, oldResultOperations, topLevelStruct2Name);

            elseif ((isstruct(field1)) && (~isstruct(field2))) || ((isstruct(field2)) && (~isstruct(field1)))
                error('One is struct the other one is not.');
            
            % Fields are cell arrays
            elseif (iscell(field1)) && (iscell(field2))
                assertCellsHaveMatchingLength(fullFieldName1, field1, fullFieldName2, field2);
                for cellIndex = 1:numel(field1)
                    cell1 = field1{cellIndex};
                    cell2 = field2{cellIndex};

                    % Update toleranceFieldPaths
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s{%d}', struct1FieldNames{fieldIndex}, cellIndex);

                    % Update list of flat struct field names
                    flatStruct1FieldNames{numel(flatStruct1FieldNames)+1} = sprintf('%s{%d}',fullFieldName1, cellIndex);
                    flatStruct2FieldNames{numel(flatStruct2FieldNames)+1} = sprintf('%s{%d}',fullFieldName2, cellIndex);

                    if (isstruct(cell1) && (isstruct(cell2))) 
                        [result, resultOperations, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames] = doRecursion(...
                            sprintf('%s{%d}',fullFieldName1, cellIndex), cell1, sprintf('%s{%d}',fullFieldName2, cellIndex), cell2, ...
                            tolerance, customTolerances, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames, ...
                            result, oldResultOperations, topLevelStruct2Name);

                    elseif (isnumeric(cell1)) && (isnumeric(cell2))
                        try
                            % See if there is a custom tolerance
                            eval(sprintf('toleranceEmployed = %f;', customTolerances(toleranceFieldPaths{numel(toleranceFieldPaths)})));
                        catch 
                            % Select default tolerance
                            toleranceEmployed = tolerance;
                        end

                        % Update the operation
                        resultOperations{numel(resultOperations)+1} = sprintf('Compare numeric arrays %s to %s using tolerance %g', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, toleranceEmployed);
                        
                        % Do the comparison
                        maxDiff = max(abs(cell1(:)-cell2(:)));
                        if (maxDiff  > toleranceEmployed)
                            result{numel(result)+1} = sprintf('Max difference between ''%s'' and ''%s'' <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, maxDiff, toleranceEmployed);
                        end

                    elseif (islogical(cell1)) && (islogical(cell2))
                        % Update the operation
                        resultOperations{numel(resultOperations)+1} = sprintf('Compare logical arrays %s to %s', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});

                        % Do the comparison
                        if (any(cell1(:)~=cell2(:)))
                            result{numel(result)+1} = sprintf('Logical arrays ''%s'' and ''%s'' do not match.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                        end

                    elseif (ischar(cell1)) && (ischar(cell2))
                        % Update the operation
                        resultOperations{numel(resultOperations)+1} = sprintf('Compare character arrays %s to %s', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});

                        % Do the comparison
                        if (~strcmp(cell1, cell2))
                            result{numel(result)+1} = sprintf('Character arrays ''%s'' and ''%s'' do not match.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                        end
                    else
                        error('Only comparing numerical, logical, char, cell and struct types\n');
                    end
                end % cellIndex

            % Fields are numeric arrays
            elseif (isnumeric(field1)) && (isnumeric(field2))

                try
                    % See if there is a custom tolerance
                    eval(sprintf('toleranceEmployed = %f;', customTolerances(toleranceFieldPaths{numel(toleranceFieldPaths)})));
                catch 
                    % Select default tolerance
                    toleranceEmployed = tolerance;
                end

                % Update the operation
                resultOperations{numel(resultOperations)+1} = sprintf('Compare numerical arrays %s to %s with tolerance %s', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, toleranceFieldPaths{numel(toleranceFieldPaths)});
            
                % Do the comparison
                maxDiff = max(abs(field1(:)-field2(:)));
                if (maxDiff  > toleranceEmployed)
                    result{numel(result)+1} = sprintf('Max difference between ''%s'' and ''%s'' <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, maxDiff, toleranceEmployed);
                end

            % Fields are logical arrays
            elseif (islogical(field1)) && (islogical(field2))
                % Update the operation
                resultOperations{numel(resultOperations)+1} = sprintf('Compare logical arrays %s to %s with tolerance %s', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, toleranceFieldPaths{numel(toleranceFieldPaths)});
            
                % Do the comparison
                if (any(field1(:)~=field2(:)))
                    result{numel(result)+1} = sprintf('Logical arrays ''%s'' and ''%s'' do not match.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                end

            % Fields are char arrays    
            elseif (ischar(field1)) && (ischar(field2))
                % Update the operation
                resultOperations{numel(resultOperations)+1} = sprintf('Compare character arrays %s to %s with tolerance %s', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, toleranceFieldPaths{numel(toleranceFieldPaths)});
            
                % Do the comparison
                if (~strcmp(field1, field2))
                    result{numel(result)+1} = sprintf('Character arrays ''%s'' and ''%s'' do not match.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                end

            % Fields are neither of the above
            else
                class(field1)
                class(field2)
                error('Only comparing numerical, logical, char, cell and struct types\n');
            end

        end
    end
end


function assertCellsHaveMatchingLength(cellArray1Name, cellArray1, cellArray2Name, cellArray2)
    if (isempty(cellArray1)) && (isempty(cellArray1))
        ;
    elseif ((isempty(cellArray1)) && (~isempty(cellArray2))) || ((isempty(cellArray2)) && (~isempty(cellArray1)))
        error('''%s'' is empty whereas ''%s'' is not. Will not compare further.', cellArray1Name, cellArray2Name);
    end

    if (numel(cellArray1) ~= numel(cellArray1))
        error('''%s'' and ''%s'' are cell arrays of different lengths (%d vs %d). Will not compare further.', cellArray1Name, cellArrayAName, numel(cellArray1), numel(cellArray2));
    end

end


function assertStructsHaveMatchingLength(struct1Name, struct1, struct2Name, struct2)

    if (isempty(struct1)) && (isempty(struct2))
        ;
    elseif ((isempty(struct1)) && (~isempty(struct2))) || ((isempty(struct2)) && (~isempty(struct1)))
        error('''%s'' is empty whereas ''%s'' is not. Will not compare further.', struct1Name, struct2Name);
    end

    % Check for non-struct inputs
    if (~isstruct(struct1)) 
        error('''%s'' is not a struct. Will not compare further.', struct1Name);
    end
    
    if (~isstruct(struct1)) 
        error('''%s'' is not a struct. Will not compare further.', struct2Name);
    end

    if (numel(struct1) ~= numel(struct2))
        error('''%s'' and ''%s'' are struct arrays of different lengths (%d vs %d). Will not compare further.', struct1Name, struct2Name, numel(struct1), numel(struct2));
    end
end

