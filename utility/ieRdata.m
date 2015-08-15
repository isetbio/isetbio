function val = ieRdata(func,rd,varargin)
% ISETBIO remote data functions
%
%   val = ieRdata(cmd,[rd],varargin)
%
% INPUTS
%  func - 'create', 'web site','file get','read image', 'dir'
%  rd   - rdata object (optional)
% 
% OUTPUT
%  val
%
% Examples:
%  rd = ieRdata('create');
%  ieRdata('web site');
%  ieRdata('dir',rd,'Stryer')
%  ieRdata('file get',[],'cText4.mat');
%  val = ieRdata('load data',rd,'cText4.mat','scene'); vcAddObject(val.scene); sceneWindow;
%  img = ieRdata('read image',rd,'birdIce'); imshow(img);
%
% BW ISETBIO Team, Copyright 2015

%%
if notDefined('func'), func = 'ls'; end
if notDefined('rd')
    % Open up the default and show the user the web-page.
    rd = rdata('base','http://scarlet.stanford.edu/validation/SCIEN');
end

%%
f = ieParamFormat(func);
switch f
    case {'create'}
        % rd = ieRdata('create')
        val = rdata('base','http://scarlet.stanford.edu/validation/SCIEN');

    case {'website'}
        rd.webSite;
        
    case {'ls','dir'}
        % dirList = ieRdata('ls',rd,'png');
        if isempty(varargin), error('Pattern required\n'); end
        val = rd.listFiles(varargin{1});
               
    case 'fileget'
        % outName = ieRdata('get',rd,fname);
        % ieRdata('fileget',[],'cText4.mat');
        %
        % One possibility.  Though maybe it should go in tempname/rdata
        if isempty(varargin), error('File string required'); end
        
        fname = varargin{1};
        destDir = fullfile(isetbioRootPath,'local');
        if ~exist(destDir,'dir')
            disp('Making local directory')
            mkdir(destDir); 
        end
        
        dest = fullfile(destDir,fname);
        val = rd.fileGet(fname,dest);
        
    case 'readimage'
        % rdata('read image', rd, fname);
        if isempty(varargin), error('remote image file name required'); end
        val = rd.readImage(varargin{1});

    case 'loaddata'
        % ieRdata('load data',rd,fname,variableName)
        %
        % Loads mat-file data, including a specific variable from a
        % matfile.
        %
        % val = ieRdata('load data',rd,'EurasianFemale_shadow.mat');
        % val = ieRdata('load data',rd,'cText4.mat','scene');
        % vcAddObject(val.scene); sceneWindow;
        if isempty(varargin), error('remote data file name required'); end
        dest = ieRdata('file get',rd,varargin{1});
        if length(varargin) == 2, val = load(dest,varargin{2});
        else val = load(dest);
        end

    otherwise
        error('Unknown ieRdata function %s\n',func);
end

end
