function val = vcConstants(con)
% Return values of physical constants used by ISET.
%
% Syntax:
%   val = vcConstants(constantName)
%
% Description:
%    We return values of various physical constants using this function, 
%    rather than a global. This will help with compilation of the Matlab
%    code. The currently stored constants are
%       {'planck', 'h', 'plancksconstant'}
%       {'q', 'electroncharge'}
%       {'c', 'speedoflight'}
%       {'j', 'joulesperkelvin', 'boltzman'}
%       {'mmPerDeg'}                         - Human retina
%
%    The code below contains examples of function usage. To access, type
%    'edit vcConstants.m' into the Command Window.
%
% Inputs:
%    con - String. The constant to retrieve.
%
% Outputs:
%    val - Numeric. The constant's value.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/10/18  jnm  Formatting

% Examples:
%{
    vcConstants('q')
    vcConstants('h')
    vcConstants('mmPerDeg')
%}

con = ieParamFormat(con);

switch lower(con)
    case {'planck', 'h', 'plancksconstant'}
        val =  6.626176e-034 ;
    case {'q', 'electroncharge'}
        val = 1.602177e-19;   % [C]
    case {'c', 'speedoflight'}
        val = 2.99792458e+8;  % Meters per second
    case {'j', 'joulesperkelvin', 'boltzman'}
        val = 1.380662e-23;  % [J/K], used in black body radiator formula
    case {'mmperdeg'}
        val = 0.3;
    otherwise
        error('Unknown physical constant');
end

return;

% Other possible constants.
%
% % electron charge:
% vcConstants.q = 1.6e-19;                       % [C]
%
% % frequency of visible lights (370nm - 730nm):
% vcConstants.nu = 3e8./((370:730)'*1e-9);       % [1/s]
%
% % Speed of light
% vcConstants.c = 2.99792458  * 10^8;            % Meters per second
%
% % Planck's Constant
% vcConstants.h =  6.6200e-034 ;
%
% % energy of photons in visible range:
% vcConstants.h_nu = 6.62e-34.*vcConstants.nu; 	 % [J]
%
% % Planck's constant:
% vcConstants.k =  1.38e-23;                     % [J/K]
%
% % Room Temperature:
% vcConstants.T = 300;                           % [K]
%
% % Permittivity of Si
% vcConstants.epsilon_s = 10.45e-13;             % [F/cm]
%
% % Permittivity of SiO2
% vcConstants.epsilon_ox = 34.5e-14;             % [F/cm]
%
% % Energy bandgap of Si:
% vcConstants.Eg = 1.12;                         % [eVolt]
