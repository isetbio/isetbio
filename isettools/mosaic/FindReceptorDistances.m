function [LDistances,MDistances,SDistances,recDistances] = FindReceptorDistances(Clist)
% [LDistances,MDistances,SDistances,recDistances] = FindReceptorDistances(Clist)
%
% For each receptor in the passed Clist, find the distance
% to the nearest L, M, and S cone in the mosaic, and to the
% nearest cone.
%
% 12/7/05   dhb     Wrote it.
% 12/11/05  dhb     Also return distance to nearest receptor.

nReceptors = size(Clist,1);
for i = 1:nReceptors
    minL = Inf; minM = Inf; minS = Inf;
    for j = 1:nReceptors
        distance = norm( Clist(i,1:2)-Clist(j,1:2) );
        if (i ~= j)
            switch (Clist(j,3))
                case 1,
                    if (distance < minL)
                        minL = distance;
                    end
                case 2,
                    if (distance < minM)
                        minM = distance;
                    end
                case 3
                    if (distance < minS)
                        minS = distance;
                    end
                otherwise
                    error('Receptor type should be 1, 2, or 3');
            end
        end
            
    end
    
    LDistances(i) = minL;
    MDistances(i) = minM;
    SDistances(i) = minS;
    recDistances(i) = min([minL minM minS]);
end