
function spikeData = rgcGetBlankSpikeDataFromFile(filename, SRG)
    % For each file
    %   Read the name, extract the cpd and contrast value
    %   Store cpd and contrast value, and train and test linTS
    % Check for missing files
    
    % format the file name if needed
    filename = strtrim(filename);
    if filename(end)~='_'
       filename = strcat(filename, '_'); 
    end
    % check the arguments : which are present?
    hasSRG = 1;
    if notDefined('SRG')
        hasSRG = 0;
    end
    
  %  for iCPD = 1:length(cpd)
   %     filenamecpd = sprintf('%scpd%d_',filename,cpd(iCPD));
    % load data into cell array
   % filename
        spikeData = cell(1,2);
       traindata = load(sprintf('%scontrastBlankTrain.mat',filename), '-mat'); 
       spikeData{1} = traindata;
       testdata = load(sprintf('%scontrastBlankTest.mat',filename), '-mat'); 
       spikeData{2} = testdata;

    %end
end