function obj = rgcCompute(obj, sensor, outersegment, varargin)
% rgcCompute: a method of @rgcLinear that computes the spiking output of the
% rgc mosaic to an arbitrary stimulus.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG


% % Compute using the superclass method
obj = rgcCompute@rgc(obj, sensor, outersegment, varargin);