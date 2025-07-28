function [fileData, theSelectedFileName] = loadFileData(varargin)
	% Parse optional input
    p = inputParser;
    p.addParameter('subDir', '', @ischar);
    p.addParameter('onlyReturnSelectedFileName', false, @islogical);
    p.parse(varargin{:});
    subdir = p.Results.subDir;

	rootDir = RGCMosaicConstructor.filepathFor.localDropboxDir();
	SLIMdataDir = 'IBIO_rgcMosaicResources/ONcenterMidgetRGCmosaics/intermediateFiles/SLIM';
	[file,path] = uigetfile({'*.mat'}, 'Select a file', fullfile(rootDir, SLIMdataDir, subdir));
	if (isempty(file))
		fileData = [];
		theSelectedFileName = '';
		return;
	end

	theSelectedFileName = fullfile(path,file);
	fileData = load(theSelectedFileName);
end
