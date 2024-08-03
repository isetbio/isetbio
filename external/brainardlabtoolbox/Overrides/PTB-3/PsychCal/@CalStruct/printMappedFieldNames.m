% Method to print the unified field names that have been mapped and their
% access property
function printMappedFieldNames(obj)

    % load all known fields
    unifiedFieldNames = keys(obj.fieldMap);
    [~,indices] = sort(lower(unifiedFieldNames));
    
    fprintf('    Available unified field names & access properties (<strong>r-o</strong> = read-only, <strong>r-w</strong> = read-write).\n');
    fprintf('    ----------------------------------------------------------------------------------------------------------------------\n');
    colsNum = 3;
    rowsNum = (round(numel(unifiedFieldNames)/colsNum)+1);
    for row = 1:rowsNum
        for col = 1:colsNum
            k = row + (col-1)*rowsNum;
            if (k <= size(unifiedFieldNames,2))
                unifiedFieldName = unifiedFieldNames{indices(k)};
                if strcmp(obj.fieldMap(unifiedFieldName).access,'read_write')
                    fprintf('%2d. %-25s (<strong>%s</strong>)        ', k, sprintf('''%s''',char(unifiedFieldName)), 'r-w');
                else
                    fprintf('%2d. %-25s (<strong>%s</strong>)        ', k, sprintf('''%s''',char(unifiedFieldName)), 'r-o');
                end
            end
        end
        fprintf('\n');
    end
    fprintf('    ----------------------------------------------------------------------------------------------------------------------\n');
end
