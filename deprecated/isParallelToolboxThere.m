function ret = isParallelToolboxThere()
% Checks wheter the parallel computing toolbox is installed
% 
%   ret = isParallelToolboxThere()
%
% (c) Stanford Synapse Team 2010

% Returns all of the toolbox names
warning('Deprecated. Use uniformed function checkToolbox instead');
vv  = ver;
ret = any(strcmp({vv.Name}, 'Parallel Computing Toolbox'));

end