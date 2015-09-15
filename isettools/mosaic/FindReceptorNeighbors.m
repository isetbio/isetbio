function [neighborIndex,nL,nM,nS] = FindReceptorNeighbors(n,coords,Clist)
% [neighborIndex,nL,nM,nS] = FindReceptorNeighbors(n,coords,Clist)
%
% Find the n receptors nearest to the passed coords, return an
% index to their rows in Clist.  Return how many are of each
% cone type.
%
% Each row of the variable coords holds a coordinate pair.
%
% 12/4/05   dhb     Wrote it.

nSpots = size(coords,1);
nReceptors = size(Clist,1);
neighborIndex = zeros(nSpots,n);
for i = 1:nSpots
    diffs = ones(nReceptors,1)*coords(i,:)-Clist(:,[2 1]);
    distances = sqrt(diffs(:,1).^2+diffs(:,2).^2);
    [nil,indices] = sort(distances);
    neighborIndex(i,:) = indices(1:n)';
end
nL = length(find(Clist(neighborIndex,3) == 1));
nM = length(find(Clist(neighborIndex,3) == 2));
nS = length(find(Clist(neighborIndex,3) == 3));
