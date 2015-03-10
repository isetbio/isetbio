function rgcCopy(h1, h2)
% for two objects of the same class h1 and h2 that are subclases of
% handle, copies all the fields of h1 in the fields of h2
% only works with properties such that (GetAccess = 'public', SetAccess = 'public')
% 
% rgcCopy(h1,h2)
% 
% ex:
% p1 = rgcParameters;
% p2 = rgcParameters;
% p1.set('name','youpi');
% rgcCopy(p1,p2);
% p2.Display('name');

if (notDefined('h1') || notDefined('h2'))
    error('You need to give two arguments');
end

if (~isequal(class(h1),class(h2)))
    error('The 2 arguments should be from the same class');
end

propNames = properties(h1);
for ii = 1:length(propNames)
    try
        h2.(propNames{ii}) = h1.(propNames{ii});
    end
end