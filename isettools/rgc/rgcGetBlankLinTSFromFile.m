function alreadyComputedBlankLinTS = rgcGetBlankLinTSFromFile(filename)

% return linTS data computed from previous simulation.
% alreadyComputedBlankLinTS{1}{layernumber} is the training data
% alreadyComputedBlankLinTS{2}{layernumber} is the test data

    % format the file name if needed
    filename = strtrim(filename);
    if filename(end)~='_'
       filename = strcat(filename, '_'); 
    end
    
    % load data into cell array
    alreadyComputedBlankLinTS = cell(1,2);
    % train data
    filename_train = strcat(filename, 'contrastBlankTrain');
    traindata = load(filename_train,'-mat');
    alreadyComputedBlankLinTS{1} = traindata.lints;
    % test data
    filename_test = strcat(filename, 'contrastBlankTest');
    testdata = load(filename_test,'-mat');
    alreadyComputedBlankLinTS{2} = testdata.lints;

end

