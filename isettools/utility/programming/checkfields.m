function bool = checkfields(s,varargin)
% Check for the presence of a field, possibly nested, in a structure
%
%  bool = checkfields(s,varargin)
%
% We often check for a field that is nested within other fields within a
% structure.   We have been doing this with a series of nested or grouped
% isfield statements. This got annoying, so I wrote this routine as a
% replacement. It extends the functionality of 'isfield'.
%
% Suppose there is a structure, pixel.OP.pd.type; You can verify that the
% sequence of nested structures is present via the call
%
%      checkfields(pixel,'OP','pd','type')
%
% A return value of 1 means the field sequence is present, 
% A return value of 0 means the sequence is absent.
%
% Copyright ImagEval Consultants, LLC, 2003.

bool  = 1;
tst   = s;

for ii = 1 : length(varargin)
    if isfield(tst, varargin{ii})
        tst = tst.(varargin{ii});
    else
        bool = 0;
        break;
    end
end

end