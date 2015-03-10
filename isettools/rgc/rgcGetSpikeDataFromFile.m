function spikeData = rgcGetSpikeDataFromFile(filename, SRG, cpd, contrast)
   % This reads spikeData from files, provided they have a regular name (it is not very robust.)
   % Only the contrast argument is actually used for now, put the other
   % arguments in the filename.
    
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
    spikeData = cell(length(contrast),2);
    for iCont=1:length(contrast)
       fprintf('contrast # %d...\n',iCont);
       traindata = load(sprintf('%scontrast%.2fTrain',filename, contrast(iCont)), '-mat'); 
       spikeData{iCont, 1} = traindata;
       testdata = load(sprintf('%scontrast%.2fTest',filename, contrast(iCont)), '-mat'); 
       spikeData{iCont, 2} = testdata;
       fprintf('...done.\n');
    end
    %end
end