%% v_rdata
%
% Test the rdata (remote data) routine
%
% This routine enables users to pull data from a remote file system on the
% web.
%
% Copyright Imageval Consulting, LLC  2015

ieInit

%% This is the base directory with SCIEN (ISET, ISETBIO, CISET) files

% We are hosting this site on scarlet.
remote.host = 'http://scarlet.stanford.edu/validation/SCIEN';

%% ls the lightfield directory for .mat files

% These are the light field data.  Because no extension is specified the
% .mat files are listed
remote.directory = fullfile('LIGHTFIELD');
rdata('cd',remote);
rdata('ls')  % Returns all the .mat files in the web page listing

%% Loading a variable from inside a matlab file

% This is an ISET scene with HDR data and a depth map
remote.directory = fullfile('LIGHTFIELD','scene');
scene = rdata('load data',remote,'benchHDR.mat','scene');

% Show the scene
vcAddObject(scene); sceneWindow;

%% Read an image file

% There are nice images here from Lubert Stryer
remote.directory = fullfile('RGB','LStryer');
img = rdata('read image',remote,'twoBirds.jpg');

% Show the image
vcNewGraphWin; imshow(img);

%% END