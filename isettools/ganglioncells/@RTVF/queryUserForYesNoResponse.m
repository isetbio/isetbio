function val = queryUserForYesNoResponse(message)
    notValidResponse = true;
    while (notValidResponse)
        message = sprintf('%s [y=YES, n=NO]: ', message);
        txt = lower(input(message, 's'));
        if (strcmp(txt, 'y')) || (strcmp(txt, 'n'))
            notValidResponse = false;
        end
    end   
    if (strcmp(txt, 'y'))
        val = true;
        fprintf('Will do.\n');
    else
        val = false;
        fprintf('Will SKIP\n');
   end
end