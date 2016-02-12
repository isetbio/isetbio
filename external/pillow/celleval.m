function v = celleval(funptr, cellarray)
%  v = celleval(funptr, cellarray)
%
%  evaluates a function 'fun' on each element of a cell array, 
%  returns a matrix the same size as cellarray
%
%  Inputs: funptr - must return a scalar
%          cellarray - cell array whose cells can be evaled by funptr
%
%  examples:  
%    cellmaxes = celleval('max', cellarray)
%    absmaxes  = celleval(@(x)max(abs(x(:))), cellarray);

csize = size(cellarray);
v = zeros(csize);
for j = 1:csize(1)
    for k = 1:csize(2)
    v(j,k) = feval(funptr, cellarray{j,k});
    end
end
