function returnString = ieKeyLicenseAnalyze
% Identifying the key-license problem
%
% ImagEval may use this routine to help customers debug problems with their
% license and key configuration.
%
% The most frequent problems are
%   * incorrectly copying and pasting key or license information
%   * invalid md5 strings.
%
% Example:
%   resturnString = ieKeyLicenseAnalyze; resturnString
%   Copy the information in resturnString and sent it in an e-mail message to
%   ImagEval.
%

%% Initialize return
returnString = 'Valid license and key';

%% Introduction for user
v = ieKeyVerify;
if isequal(v{2},1)
    fprintf('\nYour license and key are valid.\n');
    return;
else
    fprintf('\nYour license and key are not valid. ');
    fprintf('We will try to understand why.\n\n');
end

%% Check the md5 call
md5str = which('md5');
fprintf('You are using %s\n\n',md5str);
returnString = sprintf('You are using %s\n\n',md5str);

% If it is an operating-system dependent md5 (m-file), we could have a problem.
% Sometimes we don't get the right strings out of the system.
% We may also have a problem some day with 64 bit operating systems and the
% dll in Matlab, sigh.
[n,p,e] = fileparts(md5str);
if strcmp(e,'.m')
    fprintf('Operating-system dependent md5.\n\n');
    str = sprintf('Operating-system dependent md5.\n\n');
    returnString = addText(returnString,str);
    % The user has an operating system dependent md5.
    % We should test it.  The returned md5 strength length should be 32
    % characters long.
    n = length(md5('test'));
    if n ~= 32
        fprintf('md5 string length %.0d.\n\n',n);
        str = sprintf('md5 string length %.0d.\n\n',n);
        returnString = addText(returnString,str);
        return;
    end
end


%% Check that there is an ISET license
%
iLic = ieLicenseRead('iset');
if isempty(iLic)
    fprintf('** No license found.\n\n');
    str = sprintf('** No license found.\n\n');
    returnString = addText(returnString,str);

    fprintf('Now we check whether you have a key.\n');
    fprintf('You can still use ISET for a limited period of time.\n');
    fprintf('if you have a valid demo key.\n\n');
    userKey = ieKeyRead;
    if isequal(userKey{1},0)
        fprintf('** No key found.\n\n');
        str = sprintf('** No key found.\n\n');
        returnString = addText(returnString,str);
        fprintf('Please contact ImagEval for a valid license and key, or demo key.\n');
        return;
    else
        fprintf('You have a key: %s\n',userKey{1});
        str = sprintf('You have a key: %s\n',userKey{1});
        returnString = addText(returnString,str);

        keyLength = length(userKey{1});
        fprintf('It is %.0f elements long.\n',keyLength);
        str = sprintf('It is %.0f elements long.\n',keyLength);
        returnString = addText(returnString,str);

        if keyLength ~= 32
            fprintf('***The length of your key is unusual.***\n');
            str = sprintf('***The length of your key is unusual.***\n');
            returnString = addText(returnString,str);
        end
        fprintf('If ISET does not work with this key\n');
        fprintf('perhaps your demo key time has expired.\n');
        fprintf('Please contact ImagEval, and include your key in the message.\n');
    end
    return;
end

%% There is an ISET license.  So now we examine the ISET key
%
fprintf('Your ISET license is %s\n',iLic);
str = sprintf('Your ISET license is %s\n\n',iLic);
returnString = addText(returnString,str);

fprintf('We will now check for a Key.\n\n');
userKey = ieKeyRead;
if isequal(userKey{1},0)
    fprintf('** You have no key.\n\n');
    str = sprintf('** You have no key.\n\n');
    returnString = addText(returnString,str);
    fprintf('To obtain a key run ISET and follow the instructions in the\n')
    fprintf('  (Help | Get Key)\n');
    fprintf('pulldown menu in the main ISET window.\n\n');
    return;
end

fprintf('Your key is: %s\n',userKey{1});
str=sprintf('Your key is: %s\n',userKey{1});
returnString = addText(returnString,str);

fprintf('It is %.0f elements long',length(userKey{1}));
str = sprintf('It is %.0f elements long',length(userKey{1}));
returnString = addText(returnString,str);

return;



