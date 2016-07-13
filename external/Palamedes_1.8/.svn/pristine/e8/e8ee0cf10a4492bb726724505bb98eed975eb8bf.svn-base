%
%PAL_randomizeArray   Randomize array along a single dimension.
%
%syntax: [y state] = PAL_randomizeArray(x, {optional arguments})
%
%Input:
%   'x': 1, 2 or 3-D array
%
%Output:
%   'y': randomized version of input array .
%
%   'state': state of the random number generator at start of
%       randomization. Useful if randomization needs to be repeated
%       exactly. See example below.
%
%   Example: y = PAL_randomizeArray([100:2:110]) may return:
%
%       y = 
%
%           102 106 104 110 108 100
%
%   If 'x' is a 2-D or 3-D array, y = PAL_randomizeArray(x, dim) randomizes
%       'x' along dimension 'dim'. If 'dim' is not provided, randomization
%       occurs along the array's first dimension.
%
%   Example:
%
%       x = [11 12 13;
%            21 22 23;
%            31 32 33];
%
%       y = PAL_randomizeArray(x, 1) [equivalent to y = ...
%           PAL_randomizeArray(x)] randomizes the order of rows and may 
%           return:
%
%       y = 
%           
%           21  22  23
%           11  12  13
%           31  32  33
%
%       y = PAL_randomizeArray(x,2) randomizes order of columns and may 
%           return:
%
%       y = 
%           
%           12   11   13
%           22   21   23
%           32   31   33
%
%   The state of the random number generator may be set using additional
%   optional arguments. This may be useful if randomization needs to be
%   repeated exactly or to avoid unwanted exact repetition on separate
%   Matlab session (Matlab resets the random-number generator to identical 
%   states at each start-up). In order to repeat a randomization exactly,
%   record state on first randomization and set random number generator
%   state to identical state on second randomization. E.g.,:
%
%   x = [1:10];
%
%   [y state] = PAL_randomizeArray(x);
%   z = PAL_randomizeArray(x, 'state', state);
%
%   In order to set the state to a 'random' state at first call in a Matlab
%   session using the clock follow the argument 'state' by the string
%   'clock'. E.g.,:
%
%   x = [1 2 3; 4 5 6];
%   
%   [y state] = PAL_randomizeArray(x, 2, 'state', 'clock');
%
%   The above randomizes 'x' along its second dimension (i.e., it 
%   randomizes the order of the columns) after setting the random number 
%   generator to a state based on the clock and returns this state.
%
% Introduced: Palamedes version 1.0.0 (NP)
% Modified: Palamedes version 1.6.3 (see History.m)

function [x, state] = PAL_randomizeArray(x, varargin)

dim = 1;

if ~isempty(varargin)
    done = 0;
    n = 1;
    while done == 0
        if isnumeric(varargin{n})
            dim = varargin{n};
        else
            if strncmpi(varargin{n},'state',3)
                if isnumeric(varargin{n+1})
                    rand('state',varargin{n+1});
                else
                    if strncmpi(varargin{n+1},'clo',3)
                        rand('state',sum(100*clock()));
                    else
                        message = [varargin{n+1} 'is not a valid option. Ignored'];
                        warning(message);
                        done = 1;
                    end
                end
                n = n + 1;
            else
                warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});
                done = 1;
            end
        end
        if n == length(varargin)
            done = 1;
        else
            n = n + 1;
        end
    end
end

state = rand('state');

if isvector(x)
    x(1:length(x)) =  x(randperm(length(x)));
else
    enum = 1:ndims(x);
    enum(1) = dim;
    enum(dim) = 1;
    x = permute(x,enum);
    
    if ndims(x) == 2
        x(1:size(x,1),:) = x(randperm(size(x,1)),:);
    end
    if ndims(x) == 3
        x(1:size(x,1),:,:) = x(randperm(size(x,1)),:,:);
    end
    
    x = ipermute(x, enum);
            
end