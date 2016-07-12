%
%PAL_Entropy  Shannon entropy of probability density function
%
%   syntax: Entropy = PAL_Entropy(pdf,{optional argument})
%
%   Entropy = PAL_Entropy(pdf) returns the Shannon entropy (in 'nats') of 
%   the (discrete) N-D probability density function in the vector or matrix 
%   'pdf'.
%
%   Examples:
%
%   y = PAL_Entropy([.5 .5]) returns:
%
%   y = 
%       0.6931
%
%   y = log2(exp(PAL_Entropy([.5 .5]))) returns entropy in units of bits:
%
%   y = 
%       1
% 
%   y = log10(exp(PAL_Entropy([.5 .5]))) returns entropy in units of bans:
%
%   y = 
%       0.3010
%
%   Entropy = PAL_Entropy(pdf, nds), where nds is a positive integer 
%   returns an array of Entropy values calculated across the first 'nds' 
%   dimensions of 'pdf'. Passing an empty array ('[]') defaults to nds =
%   ndims(pdf);
%
%   Example:
%
%   a = [1  .5  .25;
%        0  .5  .25;
%        0  0   .25;
%        0  0   .25];
%
%   y = log2(exp(PAL_Entropy(a,1)))
%
%   y =  0     1      2
%        
%   Entropy = PAL_Entropy(pdf, nds, dm), where dm is a positive integer 
%   returns an array of Entropy values after marginalizing dimension 'dm'
%   of 'pdf'.
%
%   Example:
%
%   a = [.5 .5;
%        0  0];
%
%   y = log2(exp(PAL_Entropy(a,[],1)))
%
%   y =  1
%
%   while:
%
%   y = log2(exp(PAL_Entropy(a,[],2)))
%
%   y =  0
%        
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.1.0, 1.5.0, 1.6.0, 1.6.3 (see History.m)

function Entropy = PAL_Entropy(pdf,varargin)

nds = ndims(pdf);

if ~isempty(varargin)
    if ~isempty(varargin{1})
        nds = varargin{1};
    end
    if length(varargin) > 1
        marginalize = fliplr(sort(varargin{2}));
        for dim = 1:length(marginalize);
            pdf = sum(pdf,marginalize(dim));
        end
    end
end

Entropy = pdf.*log(pdf);
Entropy(isnan(Entropy)) = 0;          %effectively defines 0.*log(0) to equal 0.

for d = 1:nds
    Entropy = sum(Entropy,d);    
end
Entropy = -Entropy;