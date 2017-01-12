function result = RecursivelyCompareStructs(struct1Name, struct1, struct2Name, struct2, tolerance, customTolerances, graphMismatchedData, compareStringFields, result0)

    result = result0;
    
    if (isempty(struct1)) && (isempty(struct2))
        return;
    elseif (isempty(struct1)) && (~isempty(struct2))
        resultIndex = numel(result)+1;
        result{resultIndex} = sprintf('''%s'' is empty whereas ''%s'' is not. Will not compare further.', struct1Name, struct2Name);
        return;
    elseif (~isempty(struct1)) && (isempty(struct2))
        resultIndex = numel(result)+1;
        result{resultIndex} = sprintf('''%s'' is not empty whereas ''%s'' is empty. Will not compare further.', struct1Name, struct2Name);
        return;
    end
    
    % Check for non-struct inputs
    if (~isstruct(struct1)) 
        resultIndex = numel(result)+1;
        result{resultIndex} = sprintf('''%s'' is not a struct. Will not compare further.', struct1Name);
        return;
    end
    
    if (~isstruct(struct1)) 
        resultIndex = numel(result)+1;
        result{resultIndex} = sprintf('''%s'' is not a struct. Will not compare further.', struct2Name);
        return;
    end
    
    
    if (numel(struct1) ~= numel(struct2))
        resultIndex = numel(result)+1;
        result{resultIndex} = sprintf('''%s'' and ''%s'' are struct arrays of different lengths (%d vs %d). Will not compare further.', struct1Name, struct2Name, numel(struct1), numel(struct2));
        return;
    end
    
    
    % OK, inputs are good structs so lets continue with their fields
    theStruct1Name = struct1Name;
    theStruct2Name = struct2Name;
    
    for structIndex = 1:numel(struct1)   
        struct1Name = sprintf('%s(%d)', theStruct1Name, structIndex);
        struct2Name = sprintf('%s(%d)', theStruct2Name, structIndex);
        
        struct1FieldNames = sort(fieldnames(struct1(structIndex)));
        struct2FieldNames = sort(fieldnames(struct2(structIndex)));

        % Check that the two structs have same number of fields
        if numel(struct1FieldNames) ~= numel(struct2FieldNames)
            resultIndex = numel(result)+1;
            result{resultIndex} = sprintf('''%s'' has %d fields, whereas ''%s'' has %d fields. Will not compare further.', struct1Name, numel(struct1FieldNames), struct2Name, numel(struct2FieldNames));
            return;
        end
    
        for k = 1:numel(struct1FieldNames)
        
            % Check that the two structs have the same field names
            if (strcmp(struct1FieldNames{k}, struct2FieldNames{k}) == 0)
                resultIndex = numel(result)+1;
                result{resultIndex} = sprintf('''%s'' and ''%s'' have different field names: ''%s'' vs ''%s''. Will not compare further.', struct1Name, struct2Name, struct1FieldNames{k}, struct2FieldNames{k});
                return;
            end
    
            field1Name = sprintf('%s.%s', struct1Name, struct1FieldNames{k});
            field2Name = sprintf('%s.%s', struct2Name, struct2FieldNames{k});
       
            field1 = [];
            field2 = [];
            eval(sprintf('field1 = struct1(structIndex).%s;', struct1FieldNames{k}));
            eval(sprintf('field2 = struct2(structIndex).%s;', struct2FieldNames{k}));
       
           % compare structs
           if isstruct(field1)
               if isstruct(field2)
                    result = RecursivelyCompareStructs(field1Name, field1, field2Name, field2, tolerance, customTolerances, graphMismatchedData, compareStringFields, result);
               else
                    resultIndex = numel(result)+1;
                    result{resultIndex} = sprintf('''%s'' is a struct but ''%s'' is not.', field1Name, field2Name);
               end
          
               % compare strings
               elseif ischar(field1)
                   if ischar(field2)
                       if (compareStringFields)
                           if (~strcmp(field1, field2))
                                resultIndex = numel(result)+1;
                                result{resultIndex} = sprintf('''%s'' and ''%s'' are different: ''%s'' vs. ''%s''.', field1Name, field2Name, field1, field2);
                           end
                       end 
                   else
                        resultIndex = numel(result)+1;
                        result{resultIndex} = sprintf('''%s'' is a char string but ''%s'' is not.\n', field1Name, field2Name);
                   end
           
            % compare  numerics   
            elseif isnumeric(field1)
               if isnumeric(field2)
                   if (ndims(field1) ~= ndims(field2))
                       resultIndex = numel(result)+1;
                       result{resultIndex} = sprintf('''%s'' is a %d-D numeric whereas ''%s'' is a %d-D numeric.', field1Name, ndims(field1), field2Name, ndims(field2));
                   else 
                       if (any(size(field1)-size(field2)))
                            sizeField1String = sprintf((repmat('%2.0f  ', 1, numel(size(field1)))), size(field1));
                            sizeField2String = sprintf((repmat('%2.0f  ', 1, numel(size(field2)))), size(field2));
                            resultIndex = numel(result)+1;
                            result{resultIndex} = sprintf('''%s'' is a [%s] matrix whereas ''%s'' is a [%s] matrix.', field1Name, sizeField1String, field2Name, sizeField2String);
                       else
                           % equal size numerics
                           toleranceEmployed = selectToleranceToEmploy(tolerance, customTolerances, field2Name);
                           if (any(abs(field1(:)-field2(:)) > toleranceEmployed))
                                figureName = '';
                                if (graphMismatchedData)
                                    figureName = plotDataAndTheirDifference(field1, field2, field1Name, field2Name);
                                end
                                resultIndex = numel(result)+1;
                                maxDiff = max(abs(field1(:)-field2(:)));
                                if (isempty(figureName))
                                    result{resultIndex} = sprintf('Max difference between ''%s'' and ''%s'' <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>.', field1Name, field2Name, maxDiff, toleranceEmployed);
                                else
                                    result{resultIndex} = sprintf('Max difference between ''%s'' and ''%s'' <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>. See figure named: ''%s''', field1Name, field2Name, maxDiff, toleranceEmployed, figureName);
                                end
                           end
                       end
                   end
               else
                   resultIndex = numel(result)+1;
                   result{resultIndex} = sprintf('''%s'' is a numeric but ''%s'' is not.', field1Name, field2Name);
               end
           
           % compare logicals
           elseif islogical(field1)
               if islogical(field2)
                   if (ndims(field1) ~= ndims(field2))
                       resultIndex = numel(result)+1;
                       result{resultIndex} = sprintf('''%s'' is a %d-D logical whereas ''%s'' is a %d-D logical.', field1Name, ndims(field1), field2Name, ndims(field2));
                   else 
                       if (any(size(field1)-size(field2)))

                       else
                           % equal size logicals
                           if (any(field1(:) ~= field2(:)))
                               resultIndex = numel(result)+1;
                               result{resultIndex} = sprintf('There are differences between logical fields''%s'' and ''%s''.', field1Name, field2Name);
                           end
                       end
                   end
               else
                   resultIndex = numel(result)+1;
                   result{resultIndex} = sprintf('''%s'' is a logical but ''%s'' is not.', field1Name, field2Name);
               end
           
           % compare cells
           elseif iscell(field1)
               if iscell(field2)
                   if (ndims(field1) ~= ndims(field2))
                       resultIndex = numel(result)+1;
                       result{resultIndex} = sprintf('''%s'' is a %d-D cell whereas ''%s'' is a %d-D cell.', field1Name, ndims(field1), field2Name, ndims(field2));
                   else 
                       if (any(size(field1)-size(field2)))
                            sizeField1String = sprintf((repmat('%2.0f  ', 1, numel(size(field1)))), size(field1));
                            sizeField2String = sprintf((repmat('%2.0f  ', 1, numel(size(field2)))), size(field2));
                            resultIndex = numel(result)+1;
                            result{resultIndex} = sprintf('''%s'' is a [%s] matrix whereas ''%s'' is a [%s] matrix.', field1Name, sizeField1String, field2Name, sizeField2String);
                       else
                            % equal size numerics
                            result = CompareCellArrays(field1Name, field1, field2Name, field2, tolerance, customTolerances, graphMismatchedData, compareStringFields, result);
                       end
                   end
               else
                   resultIndex = numel(result)+1;
                   result{resultIndex} = sprintf('''%s'' is a cell but ''%s'' is not.', field1Name, field2Name);
               end

           % oh oh... some custom class probably.    
           else
                class(field1)
                class(field2)
                error('Do not know how to compare this class type');
            end
        end  % for k 
    end % structIndex
    
end


function result = CompareCellArrays(field1Name, field1, field2Name, field2, tolerance, customTolerances, graphMismatchedData, compareStringFields, result)

   for k = 1:numel(field1)  
       
       % Char values
       if (ischar(field1{k}))
           if (ischar(field2{k}))
               if (compareStringFields)
                   if (~strcmp(field1{k}, field2{k}))
                       resultIndex = numel(result)+1;
                       result{resultIndex} = sprintf('Corresponding cell fields have different string values: ''%s'' vs. ''%s''.', field1{k}, field2{k});
                   end
               end
           else
              resultIndex = numel(result)+1;
              result{resultIndex} = sprintf('Corresponding cell fields have different types');
           end
           
       % numeric values
       elseif (isnumeric(field1{k}))
           if (isnumeric(field2{k}))

               if (numel(field1) == 1) && (numel(field2) == 1)
                   toleranceEmployed = selectToleranceToEmploy(tolerance, customTolerances, field2Name);
                   if (abs(field1{1}-field2{1}) > toleranceEmployed)
                       resultIndex = numel(result)+1;
                       result{resultIndex} = sprintf('Corresponding cell fields have different numeric values: ''%g'' vs. ''%g''.', field1{1}, field2{1});
                   end
               else
                  subfield1 = field1{k};
                  subfield2 = field2{k};
                  if (any(size(subfield1)-size(subfield2)))
                      result{resultIndex} = sprintf('Corresponding cell subfields have different dimensionalities\n');
                  else
                      % equal size numerics
                      toleranceEmployed = selectToleranceToEmploy(tolerance, customTolerances, field2Name);
                      if (any(abs(subfield1(:)-subfield2(:)) > toleranceEmployed))
                            figureName = '';
                            if (graphMismatchedData)
                                figureName = plotDataAndTheirDifference(subfield1, subfield2, field1Name, field2Name);
                            end
                            resultIndex = numel(result)+1;
                            maxDiff = max(abs(subfield1(:)-subfield2(:)));
                            if (isempty(figureName))
                                result{resultIndex} = sprintf('Max difference between ''%s'' and ''%s'' at index %d <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>.', field1Name, field2Name, k, maxDiff, toleranceEmployed);
                            else
                                result{resultIndex} = sprintf('Max difference between ''%s'' and ''%s'' at index %d <strong>(%g)</strong> is greater than the set tolerance <strong>(%g)</strong>. See figure named: ''%s''', field1Name, field2Name, k, maxDiff, toleranceEmployed, figureName);
                            end
                       end
                  end
               end
          else
              resultIndex = numel(result)+1;
              result{resultIndex} = sprintf('Corresponding cell fields have different types');
           end
           
       % cells
       elseif (iscell(field1{k}))
           if (iscell(field2{k}))
               result = CompareCellArrays(field1Name, field1{k}, field2Name, field2{k}, tolerance, customTolerances, graphMismatchedData, compareStringFields, result);
           else
              resultIndex = numel(result)+1;
              result{resultIndex} = sprintf('Corresponding cell fields have different types');
           end
           
       % structs
       elseif isstruct(field1{k})
           if isstruct(field2{k})
                result = RecursivelyCompareStructs(field1Name, field1{k}, field2Name, field2{k}, tolerance, customTolerances, graphMismatchedData, compareStringFields, result);
           else
                resultIndex = numel(result)+1;
                result{resultIndex} = sprintf('''%s'' is a struct but ''%s'' is not.', field1Name, field2Name);
           end
        
       % oh oh... some custom class probably.
       else
          class(field1{k})
          class(field2{k})
          error(2,'structsAreSimilar.CompareCellArrays. non-char, non-numeric comparison not implemented yet\n');
       end
   end
end


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
 

% Method to select which tolerance to employ by checking if the fieldName exists in customTolerances
function toleranceEmployed = selectToleranceToEmploy(globalTolerance, customTolerances, fieldName)
    toleranceEmployed = globalTolerance;
    
    fieldName = strrep(fieldName, 'validationData.', '');
    drilledDownToleranceStuct = customTolerances;
    drilledFieldName = fieldName;
    
    fieldPath = strsplit(drilledFieldName, '.');
    for k = 1:numel(fieldPath)
        if (isfield(drilledDownToleranceStuct, fieldPath{k}))
            drilledDownToleranceStuct = drilledDownToleranceStuct.(fieldPath{k});
            drilledFieldName = strrep(drilledFieldName, sprintf('%s.', fieldPath{k}), '');
            if (k ==numel(fieldPath))
                toleranceEmployed = drilledDownToleranceStuct;
            end
        end
    end

    %fprintf('Tolerance for field %s: %g\n', fieldName, toleranceEmployed);
end
