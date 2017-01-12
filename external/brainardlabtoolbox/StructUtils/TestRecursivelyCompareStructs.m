function TestRecursivelyCompareStructs()

    % Generate a struct
    s1.aStructArray = [struct('f1', 1, 'f2', 'as'), struct('f1', -1, 'f2', 'as')];
    s1.aCellArray = {1, 'o', [1 2 3]};
    s1.aLogicalArray = [true true false];
    s1.aNumericalArray =  [ -1 0 1000 inf eps];
    s1.aString = 'just a test';
    s1.aFloat = 12.34;
    s1.aNestedCellArray = {'n1', struct('f1', 1, 'f2', true, 'f3', struct('sub1', 2)), [1 2 3]};
    
    % Make the second struct identical to the first one
    s2 = s1;
    
    % Unit tests
    % Difference in a struct field within an array of structs
    %s2.aStructArray(2) = struct('f1', 1, 'f2', 'bs');
    
    % Cell array with different field types
    %s2.aCellArray = {1, 2, [1 2 3]};
    
    % Difference in a logical array
    %s2.aLogicalArray = [true true true];
    
    % Difference in a numerical array
    %s2.aNumericalArray(3) = s2.aNumericalArray(3) + 1.1e-12;
    
    % Differences in a nested cell array
    %s2.aNestedCellArray{1} = 'n2';
    %s2.aNestedCellArray{2}.f1 = 1.1;
    %s2.aNestedCellArray{2}.f3.sub1 = 2.01;
    
    compareStringFields = true;
    graphMismatchedData = false;
    customTolerances = [];
    tolerance = 1e-12;
    
    
    % Call the function
    result = RecursivelyCompareStructs('s1', s1, 's2', s2, tolerance, customTolerances, graphMismatchedData, compareStringFields, []);
    
    % Print the results
    if (isempty(result))
        fprintf('\n<strong>Structs comparison results: structs are identical.</strong>\n');
    else
        fprintf('\n<strong>Structs comparison results: structs are NOT identical. </strong>\n');
        for k = 1:numel(result)
            fprintf('\t[%d]. %s\n', k, result{k});
        end
    end
end

