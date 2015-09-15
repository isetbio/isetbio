%% t_codeRemoteData
%
%  Tutorial
%
% BW ISETBIO Team, Copyright 2015

%
ieInit

% In ISETBIO we create a Matlab object to manage remote data.  The
% interface is
rd = ieRdata('create');

% This object has pointers to the URL where the data are stored as well as
% a list of files and directories at that site.

% Here are the slots in the object
rd

% The data are summarized as a list of directories and for each directory
% there is a list of files

% Here is a dump of the directories and the files
ieRdata('files print',rd)

% These were created on the remote site using the utility in the
% RemoteDataToolbox  rdSiteTOC.m
%
% This utility is run on black, usually, and it walks the directory tree to
% find all the directories and files.  It stores them in a Table of
% Contents in the directory it is run from.
%
% The TOC file is a JSON file.  That is the file we read to assign the
% slots in the rd object.
%

% To view the web-site directly you can use
ieRdata('web site',rd);

%% To download a a directory listing from the remote site use 

val =  ieRdata('ls',rd,'Stryer');
disp(val)

%% To download a file use
localFile = ieRdata('file get',rd,'cText4.mat')

%% You can load a single variable from a remote Matlab file
val = ieRdata('load data',rd,'cText4.mat','scene');
disp(val)
ieAddObject(val.scene); sceneWindow;

%% You can read an image (png, jpg, or other file types) as
img = ieRdata('read image',rd,'birdIce');
vcNewGraphWin; imshow(img);

%%
