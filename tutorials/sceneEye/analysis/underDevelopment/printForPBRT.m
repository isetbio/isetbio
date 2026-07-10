function printForPBRT(v)

for ii = 1:length(v)
    
    if(mod(v(ii),1) == 0)
        fprintf('%i, ',v(ii));
    else
        fprintf('%0.7f, ',v(ii));
    end
    
    if(mod(ii,10) == 0 || ii == length(v))
        fprintf('\n');
    end
end

end
