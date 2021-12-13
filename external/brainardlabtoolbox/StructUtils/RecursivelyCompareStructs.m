function result = RecursivelyCompareStructs(struct1Name, struct1, struct2Name, struct2, varargin)
% result = RecursivelyCompareStructsTests(struct1Name, struct1, struct2Name, struct2, varargin)
%
% Method to compare nested structs with arbitrary internal organization.
% Key/value pairs
% 'defaultTolerance', a numeric value applied to all differentials
% 'customTolerance', a containers.Map() with tolerance values for specific fields (see RecursivelyCompareStructsTests).
% 'graphMismatchedData', boolean, specifying wether to plot any mismatched (numerical only) arrays.
%
% Simple usage:
% s1 = struct('a', 1);
% s2 = struct('a', 2);
% result = RecursivelyCompareStructs('s1', s1, 's2', s2)
%
% For more elabored usage examples see RecursivelyCompareStructsTests.
%
% 02/07/17  NPC  Wrote it.
% 10/05/17  NPC  If a field is not found in one struct do not throw an
%                error. Print a message in red and skip to the next field.
%% Parse input
p = inputParser;
p.addParameter('defaultTolerance', 0, @isnumeric);
p.addParameter('customTolerances', containers.Map());
p.addParameter('graphMismatchedData', true);
p.addParameter('failOnMissingFields', false);
p.addParameter('verbosityLevel', 1);
p.parse(varargin{:});
defaultTolerance = p.Results.defaultTolerance;
customTolerances = p.Results.customTolerances;
graphMismatchedData = p.Results.graphMismatchedData;
failOnMissingFields = p.Results.failOnMissingFields;
verbosityLevel = p.Results.verbosityLevel;
% Go
[result, resultOperations, ~, ~, ~] = doRecursion(...
    struct1Name, struct1, struct2Name, struct2, ...
    defaultTolerance, customTolerances, {}, {}, {}, ...
    [], [], struct2Name, graphMismatchedData, failOnMissingFields, verbosityLevel);

% Print the results
if (verbosityLevel > 0)
    if (isempty(result))
        %fprintf('\n<strong>Structs comparison results: structs are identical.</strong>\n');
    else
        fprintf('\n<strong>Structs comparison results: structs are NOT identical. </strong>\n');
        for k = 1:numel(result)
            fprintf('\t[%d]. %s\n', k, result{k});
        end
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
    oldResult, oldResultOperations, topLevelStruct2Name, graphMismatchedData, failOnMissingFields, verbosityLevel)

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

        if (numel(struct1FieldNames) ~= numel(struct2FieldNames))
            fprintf(2,'** structs ''%s'' and ''%s'' have different number of fields: %d vs %d\n', struct1Name, struct2Name, numel(struct1FieldNames), numel(struct2FieldNames))
            if (failOnMissingFields)
                result{numel(result)+1} = sprintf('** structs ''%s'' and ''%s'' have different number of fields: %d vs %d\n', struct1Name, struct2Name, numel(struct1FieldNames), numel(struct2FieldNames));
                continue;
            end
        end
        
        for fieldIndex = 1:numel(struct1FieldNames)
            theFieldName = struct1FieldNames{fieldIndex};
            fullFieldName1 = sprintf('%s.%s', struct1FullName,theFieldName);
            if (~isfield(struct2, theFieldName))
                if (failOnMissingFields)
                    result{numel(result)+1} = sprintf('*** Field ''%ss'' not found in struct ''%s''. Skipping this field.\n', theFieldName, struct2FullName);
                    continue;
                else
                    fprintf(2, '**** Field ''%ss'' not found in struct ''%s''. Skipping this field.\n', theFieldName, struct2FullName);
                    continue;
                end
            end
            
            fullFieldName2 = sprintf('%s.%s', struct2FullName,theFieldName);
            field1 = [];
            field2 = [];
            eval(sprintf('field1 = struct1(structIndex).%s;', theFieldName));
            eval(sprintf('field2 = struct2(structIndex).%s;', theFieldName));

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
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s', theFieldName);
                else
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s.%s', structToleranceName, theFieldName);
                end
            end

            % Fields are structs
            if (isstruct(field1)) && (isstruct(field2))
                % Go into recursion
                [result, resultOperations, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames] = doRecursion(...
                    fullFieldName1, field1, fullFieldName2, field2, ...
                    tolerance, customTolerances, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames, ...
                    result, oldResultOperations, topLevelStruct2Name, graphMismatchedData);

            elseif ((isstruct(field1)) && (~isstruct(field2))) || ((isstruct(field2)) && (~isstruct(field1)))
                error('One is struct the other one is not.');
            
            % Fields are cell arrays
            elseif (iscell(field1)) && (iscell(field2))
                assertCellsHaveMatchingLength(fullFieldName1, field1, fullFieldName2, field2);
                for cellIndex = 1:numel(field1)
                    cell1 = field1{cellIndex};
                    cell2 = field2{cellIndex};

                    % Update toleranceFieldPaths
                    toleranceFieldPaths{numel(toleranceFieldPaths)+1} = sprintf('%s{%d}', theFieldName, cellIndex);

                    % Update list of flat struct field names
                    flatStruct1FieldNames{numel(flatStruct1FieldNames)+1} = sprintf('%s{%d}',fullFieldName1, cellIndex);
                    flatStruct2FieldNames{numel(flatStruct2FieldNames)+1} = sprintf('%s{%d}',fullFieldName2, cellIndex);

                    if (isstruct(cell1) && (isstruct(cell2))) 
                        [result, resultOperations, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames] = doRecursion(...
                            sprintf('%s{%d}',fullFieldName1, cellIndex), cell1, sprintf('%s{%d}',fullFieldName2, cellIndex), cell2, ...
                            tolerance, customTolerances, toleranceFieldPaths, flatStruct1FieldNames, flatStruct2FieldNames, ...
                            result, oldResultOperations, topLevelStruct2Name, graphMismatchedData);

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
                            if (graphMismatchedData)
                                figureName = plotDataAndTheirDifference(cell1, cell2, flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                            end
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
                    if (graphMismatchedData)
                       figureName = plotDataAndTheirDifference(field1, field2, flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)});
                    end
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

            % Fields  do not match
            elseif (~strcmp(class(field1), class(field2)))
                result{numel(result)+1} = sprintf('fields ''%s'' and ''%s'' have non-matching types: ''%s'' vs ''%s''.', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, flatStruct2FieldNames{numel(flatStruct2FieldNames)}, class(field1), class(field2));
            
            % Fields are neither of the above
            else
                    fprintf(2,'Class of field ''%s'': ''%s''\n', flatStruct1FieldNames{numel(flatStruct1FieldNames)}, class(field1));
                    fprintf(2,'Class of field ''%s'': ''%s''\n', flatStruct1FieldNames{numel(flatStruct2FieldNames)}, class(field2));
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

% Plotting subroutines

% Method to plot mistmatched data
function figureName = plotDataAndTheirDifference(field1, field2, field1Name, field2Name)
  
    h = figure();
    figureName = sprintf('''%s'' vs. ''%s''', field1Name, field2Name);
    set(h, 'Name', figureName);
    clf;
    
    if (numel(find(size(field1)>1)) == 1)
        % essentially a 1D vector
        field1 = field1(:);
        field2 = field2(:);
        
        diff = field1 - field2;
        minAll = min([min(field1(:)) min(field2(:))]);
        maxAll = max([max(field1(:)) max(field2(:))]);
        maxDiff = max(abs(diff(:)));
        
        set(h, 'Position', [100 100 1400 380]);
        subplot(1,3,1);
        plot(1:numel(field1), field1, 'r-');
        set(gca, 'YLim', [minAll maxAll]);
        title(sprintf('''%s''', field1Name));

        subplot(1,3,2);
        plot(1:numel(field2), field2, 'b-');
        set(gca, 'YLim', [minAll maxAll]);
        title(sprintf('''%s''', field2Name));

        subplot(1,3,3);
        plot(1:numel(field1), field1-field2, 'k-');
        set(gca, 'YLim', maxDiff * [-1 1]);
        title(sprintf('''%s'' - \n''%s''', field1Name, field2Name));
        
    elseif (ndims(field1) == 2) && (all(size(field1) > 1))
        diff = field1 - field2;
        minAll = min([min(field1(:)) min(field2(:))]);
        maxAll = max([max(field1(:)) max(field2(:))]);
        maxDiff = max(abs(diff(:)));
        
        set(h, 'Position', [100 100 1400 380]);
        subplot(1,3,1);
        imagesc(field1);
        colorbar
        set(gca, 'CLim', [minAll maxAll]);
        title(sprintf('''%s''', field1Name));

        subplot(1,3,2);
        imagesc(field2);
        set(gca, 'CLim', [minAll maxAll]);
        colorbar
        title(sprintf('''%s''', field2Name));

        subplot(1,3,3);
        imagesc(diff);
        set(gca, 'CLim', [-maxDiff maxDiff]);
        colorbar
        title(sprintf('''%s'' - \n''%s''', field1Name, field2Name));
        colormap(gray(256));
        
    elseif (ndims(field1) == 1) || ((ndims(field1)==2) && (any(size(field1)==1)))  
        diff = field1 - field2;
        minAll = min([min(field1(:)) min(field2(:))]);
        maxAll = max([max(field1(:)) max(field2(:))]);
        maxDiff = max(abs(diff(:)));
        delta = (maxAll-minAll)/10;
        
        set(h, 'Position', [100 100 1400 380]);
        subplot(1,3,1);
        plot(field1, 'bs-', 'MarkerFaceColor', [0.8 0.8 1]);
        set(gca, 'YLim', [minAll-delta maxAll+delta]);
        title(sprintf('''%s''', field1Name));
        
        subplot(1,3,2);
        plot(field2, 'bs-',  'MarkerFaceColor', [0.8 0.8 1]);
        set(gca, 'YLim', [minAll-delta maxAll+delta]);
        title(sprintf('''%s''', field2Name));
        
        subplot(1,3,3);
        plot(field1 - field2, 'bs-',  'MarkerFaceColor', [0.8 0.8 1]);
        set(gca, 'YLim', double(maxDiff)*[-1.1 1.1]);
        title(sprintf('''%s'' - \n''%s''', field1Name, field2Name));

    elseif (ndims(field1) == 3)
        diff = field1 - field2;
        minAll = min([min(field1(:)) min(field2(:))]);
        maxAll = max([max(field1(:)) max(field2(:))]);
        maxDiff = max(abs(diff(:)));
        
        [d1, d2, d3] = size(field1);
        if (d3 > 50)
            figSize = [2350 810];
            showTicks = false;
            wMargin = 0.005;
        else
            figSize = [1800 800];
            showTicks = true;
            wMargin = 0.01;
        end
        set(h, 'Position', [100 100  figSize(1)  figSize(2)]);
        
        halfD = round(d3/2);
        posVectors = getSubPlotPosVectors(...
            'rowsNum', 6, 'colsNum', halfD, ...
            'widthMargin', wMargin, 'heightMargin', 0.01, ...
            'leftMargin', 0.02', ...
            'bottomMargin', 0.01, 'topMargin', 0.01);

        for k = 1:halfD
           row = 1;
           col = k;
           subplot('Position', posVectors(row,col).v);
           
            if (~isreal(field1))
               warndlg('Displaying the ABS(field2)', 'Complex Variable !!');
               field2 = abs(field1);
           end
           imagesc(squeeze(field1(:,:,k)));
           if (k == 1)
               ylabel(field1Name, 'Color', [1 0 0]);
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end
           set(gca, 'CLim', [minAll maxAll]);
           title(sprintf('i=%d', k), 'Color', [1 0 0]);
           %colorbar('horiz');
           
           row = 2;
           col = k;
           subplot('Position', posVectors(row,col).v);
           
           if (~isreal(field2))
               warndlg('Displaying the ABS(field2)', 'Complex Variable !!');
               field2 = abs(field2);
           end
           imagesc(squeeze(field2(:,:,k)));
           if (k == 1)
               ylabel(field2Name, 'Color', [0 0 1]);
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end
           set(gca, 'CLim', [minAll maxAll]);
           title(sprintf('i=%d', k), 'Color', [0 0 1]);
           %colorbar('horiz');
           
           row = 3;
           col = k;
           subplot('Position', posVectors(row,col).v);
           if (~isreal(diff))
               warndlg('Displaying the ABS(diff)', 'Complex Variable !!');
               diff = abs(diff);
           end
           imagesc(squeeze(diff(:,:,k)));
           if (k == 1)
               ylabel(sprintf('groundTruthData\n - validationData'));
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end

           set(gca, 'CLim', [-maxDiff maxDiff]);
           title(sprintf('i=%d', k));
           %colorbar('horiz');
        end
        
        for k = halfD+1:d3
           row = 4;
           col = k-halfD;
           subplot('Position', posVectors(row,col).v);
           imagesc(squeeze(field1(:,:,k)));
           if (k == halfD+1)
               ylabel(field1Name, 'Color', [1 0 0]);
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end
           set(gca, 'CLim', [minAll maxAll]);
           title(sprintf('i=%d', k), 'Color', [1 0 0]);
           %colorbar('horiz');
           
           row = 5;
           col = k-halfD;
           subplot('Position', posVectors(row,col).v);
           imagesc(squeeze(field2(:,:,k)));
           if (k == halfD+1)
               ylabel(field2Name, 'Color', [0 0 1]);
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end
           set(gca, 'CLim', [minAll maxAll]);
           title(sprintf('i=%d', k), 'Color', [0 0 1]);
           %colorbar('horiz');
           
           row = 6;
           col = k-halfD;
           subplot('Position', posVectors(row,col).v);
           imagesc(squeeze(diff(:,:,k)));
           if (k == halfD+1)
               ylabel(sprintf('groundTruthData\n - validationData'));
           end
           axis 'image'
           if (~showTicks)
              set(gca, 'XTick', [], 'YTick', []);
           end
           set(gca, 'CLim', [-maxDiff maxDiff]);
           title(sprintf('i=%d', k));
        end
        colormap(gray(512));
    end
end

function posVectors = getSubPlotPosVectors(varargin)

    self.rowsNum        = 2;
    self.colsNum        = 2;
    self.widthMargin    = 0.01;
    self.heightMargin   = 0.01;
    self.leftMargin     = 0.06;
    self.bottomMargin   = 0.08;
    self.topMargin      = 0.08;
    
    % parse inputs
    parser = inputParser;
    parser.addParamValue('rowsNum',         self.rowsNum);
    parser.addParamValue('colsNum',         self.colsNum);
    parser.addParamValue('widthMargin',     self.widthMargin);
    parser.addParamValue('heightMargin',    self.heightMargin);
    parser.addParamValue('leftMargin',      self.leftMargin);
    parser.addParamValue('bottomMargin',    self.bottomMargin); 
    parser.addParamValue('topMargin',       self.topMargin);

    
    % Execute the parser to make sure input is good
    parser.parse(varargin{:});
    pNames = fieldnames(parser.Results);
    for k = 1:length(pNames)
       self.(pNames{k}) = parser.Results.(pNames{k}); 
    end
    
    plotWidth  = ((1.0-self.leftMargin) - self.widthMargin*(self.colsNum-1) - 0.01)/self.colsNum;
    plotHeight = ((1.0-self.bottomMargin-self.topMargin) - self.heightMargin*(self.rowsNum-1) - 0.01)/self.rowsNum;
    
    for row = 1:self.rowsNum
        yo = 0.99 - self.topMargin - (row)*(plotHeight+self.heightMargin) + self.heightMargin;
        for col = 1:self.colsNum
            xo = self.leftMargin + (col-1)*(plotWidth+self.widthMargin);
            posVectors(row,col).v = [xo yo plotWidth plotHeight];
        end
    end
end