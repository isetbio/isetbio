function b = ieReadBoolean(question)
% Get a yes or no answer. If only the world were this simple.
%
% Syntax:
%   b = ieReadBoolean(question)
%
% Description:
%    A Text dialogue box which displays the provided (or default) question
%    to the user with the optional buttons 'Yes', 'No', and 'Cancel'
%    representing the Boolean Options True (Yes), False (No), and Empty
%    (Cancel).
%
%    Please note that proposing through this dialogue box is ill-advised.
%
%    There are examples contained in the code below. To access, type 'edit
%    ieReadBoolean.m' into the Command Window.
%
% Inputs:
%    question - (Optional) String. The dialog to display to the user in the
%               appearing text box. Default is 'Yes or No?'
%
% Outputs:
%    b        - Boolean. The boolean answer in numerical Format. Options
%               are [], 0, 1.
%
% Optional key/value pairs:
%    None.
%

% Examples:
%{
    % ETTBSkip - Windows require user input.
    b = ieReadBoolean('Will you marry me?')
    b = ieReadBoolean('Have you seen Firefly?')
    b = ieReadBoolean('Have you done your homework?')
%}

if notDefined('question'), question = 'Yes or No?'; end
b = [];

ButtonName = questdlg(question, 'Options', 'Yes', 'No', 'Cancel', 'Yes');
if isempty(ButtonName), return; end

switch ButtonName
    case 'Yes'
        b = 1;
    case 'No'
        b = 0;
    case 'Cancel'
        b = [];
        return;
end

end