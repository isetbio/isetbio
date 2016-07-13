%
% PAL_sqDistanceYfuncX
%   
%   Internal function
%
% Introduced: Palamedes version 1.0.0 (FK)

function sqDistance=PAL_sqDistanceYfuncX(X,Y,func,params, varargin)

if length(params)==0
sqDistance=(Y-func(X,varargin{:})).^2;
end

if length(params)==1
sqDistance=(Y-func(X,params(1),varargin{:})).^2;
end

if length(params)==2
sqDistance=(Y-func(X,params(1),params(2),varargin{:})).^2;
end

if length(params)==3
sqDistance=(Y-func(X,params(1),params(2),params(3),varargin{:})).^2;
end