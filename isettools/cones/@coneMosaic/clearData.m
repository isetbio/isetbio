function obj = clearData(obj, varargin)
%CLEARDATA  Clear computed data
%    obj = CLEARDATA(obj, varargin)
%
%    Clear out data from object that is computed.  Useful to prevent use of
%    stale data.

% HJ ISETBIO Team 2016

    obj.absorptions = [];
    obj.current = [];
end
