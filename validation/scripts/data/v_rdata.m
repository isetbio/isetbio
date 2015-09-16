function varargout = v_rdata(varargin)
%
% Test the rdata (remote data) routine
%
% This routine enables users to pull data from a remote file system on the
% web.
%
% Copyright Imageval Consulting, LLC  2015
%
% 7/15/15  dhb  Brought this into the UnitTestToolbox world.

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);

end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Initialize ISETBIO
close all; ieInit;

%% Some informative text
UnitTest.validationRecord('SIMPLE_MESSAGE', 'Validate rdata function.');

%% This is the base directory with SCIEN (ISET, ISETBIO, CISET) files

rd = ieRdata('create'); % Test that it opens

% ieRdata('web site');    % Open the web page.

val =  ieRdata('dir',rd,'Stryer');
disp(val)

ieRdata('file get',[],'cText4.mat');

val = ieRdata('load data',rd,'cText4.mat','scene');
disp(val)
% if (runTimeParams.generatePlots)
%     vcAddObject(val.scene); sceneWindow;
% end

img = ieRdata('read image',rd,'birdIce');
if (runTimeParams.generatePlots)
    vcNewGraphWin; imshow(img);
end

end