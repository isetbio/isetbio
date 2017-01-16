function clearData(obj, varargin)
%CLEARDATA  Clear computed data
%   CLEARDATA(obj, varargin)
%
%   Clear out data from object that is computed.  Useful to prevent use of
%   stale data.
%
%   Inputs:
%   obj - a coneMosaic object

% HJ ISETBIO Team 2016

    obj.absorptions = [];
    obj.current = [];
end
