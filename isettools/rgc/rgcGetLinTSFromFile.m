
function alreadyComputedLinTS = rgcGetLinTSFromFile(filename, SRG, cpd, contrast)
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
    hasCpd = 1;
    if notDefined('cpd')
        hasCpd = 0;
    end
    hasContrast = 1;
    if notDefined('contrast')
        hasContrast = 0;
    end
    
  %  for iCPD = 1:length(cpd)
   %     filenamecpd = sprintf('%scpd%d_',filename,cpd(iCPD));
    % load data into cell array
   % filename
    alreadyComputedLinTS = cell(length(contrast),2);
    for iCont=1:length(contrast)

       traindata = load(sprintf('%scontrast%.2fTrain',filename, contrast(iCont)), '-mat'); 
       alreadyComputedLinTS{iCont, 1} = traindata.lints;
       testdata = load(sprintf('%scontrast%.2fTest',filename, contrast(iCont)), '-mat'); 
       alreadyComputedLinTS{iCont, 2} = testdata.lints;
    end
    %end
end