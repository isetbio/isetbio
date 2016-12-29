function status = isemptyorstruct(value)
% status = isemptyorstruct(value)
% 
% Returns true if the input is empty or a structure

if (isempty(value))
    status = true;
    return;
end

if (isstruct(value))
    status = true;
    return;
end