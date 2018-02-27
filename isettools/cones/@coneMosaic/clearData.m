function clearData(obj, varargin)
% Clear computed data
%
% Syntax:
%   CLEARDATA(obj, varargin)
%
% Description:
%    Clear out data from object that is computed.  Useful to prevent use of
%    stale data.
%
% Inputs:
%    obj - a coneMosaic object
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/22/18  jnm  Formatting

    obj.absorptions = [];
    obj.current = [];
end
