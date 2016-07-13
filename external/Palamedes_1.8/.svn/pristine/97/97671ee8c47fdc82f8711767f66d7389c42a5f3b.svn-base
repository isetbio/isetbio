%
% PAL_minimize  Nelder-Mead simplex function minimization
%
%   syntax: [ x fval exitflag output] = PAL_minimize( fun, x0, options, 
%       varargin )
%
%   Internal function
%
%   Performs Nelder-Mead simplex minimization of function.
%   input-output structure of this function is similar to Matlab's 
%       fminsearch function. Fields Display, TolX, TolFun, MaxIter, 
%       MaxFunEvals in options structure similar to fminsearch (though they 
%       behave somewhat different and differ in their default values). 
%
%   Notable difference between PAL_minimize and fminsearch:
%       -default tolerance values are 1e-6 (1e-4 in fminsearch)
%       -x-guesses need to differ from 0 by minimum amount before initial
%           simplex side is made proportional to value. This avoids value 
%           from getting stuck near 0 when initial guess is (very) near 
%           zero.
%
%   Standard Nelder-Mead algorithm (e.g., Nelder, J.A. & Mead, R. (1965). A
%   simplex method for function minimization. Computer Journal, 7,
%   308-313).
%
% Introduced: Palamedes version 1.4.0 (NP)
% Modified: Palamedes version 1.6.2, 1.6.3 (see History.m)

function [ x, fval, exitflag, output] = PAL_minimize( fun, x0, options, varargin )

    
%Return options structure with default and display information:
if nargin<=2 && nargout <= 1 && strcmpi(fun,'options')
    options = struct('Display','off','MaxIter','400*numberOfVariables',...
    'MaxFunEvals','400*numberOfVariables','TolX',1e-6,'TolFun',1e-6);
    if nargin == 2 && strcmpi(x0,'help')
        fprintf('\nNote that the field names are case sensitive!\n\n');
        disp(options);
        fprintf('\nDisplay:\n\t''off'': no output.\n\t''notify'': message is ');
        fprintf('displayed only in case search does not converge.\n\t''final'': ');
        fprintf('message is displayed after search is terminated.\n\t''iter'': ');
        fprintf('message is displayed after each iteration.\n\nMaxIter:\n\t');
        fprintf('maximum number of search iterations to perform. Must either be the\n\t');
        fprintf('exact default string or an integer.\n\nMaxFunEvals:\n\tmaximum number ');
        fprintf('of function evaluations to perform. Must either be the\n\texact default ');
        fprintf('string or an integer.\n\nTolX:\n\ttolerance on parameter values.\n\n');
        fprintf('TolFun:\n\ttolerance on function values.\n\n');
    end
    x = options;
    return
end

n = length(x0);

tolX = 1e-6;
tolFun = 1e-6;
maxIter = 400*n;
maxFunEvals = 400*n;
display = 0;

if ~isempty(options)
    if isfield(options,'TolX')
        tolX = options.TolX;
    end
    if isfield(options,'TolFun')
        tolFun = options.TolFun;
    end
    if isfield(options,'MaxIter')
        maxIter = options.MaxIter;
    end
    if isfield(options,'MaxFunEvals')
        maxFunEvals = options.MaxFunEvals;
    end    
    if isfield(options,'Display')
        switch lower(options.Display)
            case 'off'
                display = 0;
            case 'notify'
                display = 1;
            case 'iter'
                display = 2;
            case 'final'
                display = 3;
        end
    end            
end

if ischar(maxFunEvals)
    if ~strcmpi(maxFunEvals,'400*numberOfVariables');
        warningMessage = strcat(maxFunEvals,' is not a valid option for maxFunEvals. Use a scalar value or (the string) ''400*numberofvariables''.');
        if strcmpi(maxFunEvals,'200*numberOfVariables');
            warningMessage = strcat(warningMessage,' You might be using default values for Matlab''s fminsearch. Use: options = PAL_minimize(''options''); instead of: options = optimset(''fminsearch'');');
        end
        warning('PALAMEDES:invalidOption', warningMessage);
    end
    maxFunEvals = 400*n;
end
if ischar(maxIter)      
    if ~strcmpi(maxIter,'400*numberOfVariables');
        warningMessage = strcat(maxIter,' is not a valid option for maxIter. Use a scalar value or (the string) ''400*numberofvariables''.');
        if strcmpi(maxIter,'200*numberOfVariables');
            warningMessage = strcat(warningMessage,' You might be using default values for Matlab''s fminsearch. Use: options = PAL_minimize(''options''); instead of: options = optimset(''fminsearch'');');
        end
        warning('PALAMEDES:invalidOption', warningMessage);
    end
    maxIter = 400*n;
end

V = zeros(n+1,n);   %Simplex verteces
f = zeros(n+1,1);   %function values at verteces
V(1,:) = x0;        

iterations = 0;
evals = 0;
exitflag = logical(false);

%create initial simplex and evaluate function at verteces
%credit L. Pfeffer at Stanford
delta = 0.05;   %determines initial simplex size
    
for j = 1:n
    V(j+1,:) = x0;
    if abs(V(j+1,j)) < delta*.00025 %this (as opposed to '== 0') avoids trouble 
                                %in cases where y(j) differs by a tiny 
                                %amount from 0 (e.g., due to rounding error)
        V(j+1,j) = .00025;
    else
        V(j+1,j) = (1 + delta)*V(j+1,j);    
    end    
end
for j = 1:n+1
    f(j) = fun(V(j,:),varargin{:});
end

if display == 2
    fprintf('\titer:\tevals:\tfeval:\tparams:\n');
    fprintf('%6d\t%6d\t%9.3f\t',iterations,evals,f(1));     
    for j = 1:n            
        fprintf('%9.3f\t',V(1,j));    %user-supplied guesses
    end
    fprintf('user guess\n');
end

iterations = 1;
evals = n+1;

[f, I] = sortrows(f);
V = V(I,:);

if display == 2
    fprintf('%6d\t%6d\t%9.3f\t',iterations,evals,f(1));
    for j = 1:n            
        fprintf('%9.3f\t',V(1,j));
    end
    fprintf('initial simplex\n');
end

%Nelder-Mead search params
paramsNM = [1 2 .5 .5]; %reflection expansion contraction shrink (standard values)

%tolX criterion: longest distance across all but the best vertex to the 
%   best vertex as projected unto coordinate axis in parameter space.
%tolFun criterion: greatest difference between lowest function value and
%   others in simplex.

%Nelder-Mead iteration:

while (max(max(abs(V(2:n+1,:)-V(ones(1,n),:)))) > tolX || max(abs(f(2:n+1)-f(1))) > tolFun) && iterations < maxIter && evals < maxFunEvals        

    shrink = 0;
    centroid = mean(V(1:n,:),1);

    xr = centroid + paramsNM(1)*(centroid - V(n+1,:)); %reflection
    fr = fun(xr,varargin{:});
    evals = evals + 1;
    
    if fr < f(n) && fr >= f(1)
        Vnew = xr;
        fnew = fr;
        how = 'reflection';
    else
        if fr < f(1)
            xe = centroid + paramsNM(2)*(xr - centroid);  %expansion        
            fe = fun(xe,varargin{:});
            evals = evals + 1;
            if fe < fr
                Vnew = xe;
                fnew = fe;
                how = 'expansion';
            else
                Vnew = xr;
                fnew = fr;
                how = 'reflection';
            end
        else            
            if fr < f(n+1) && fr >= f(n)
                xoc = centroid + paramsNM(3)*(xr - centroid); %outside contraction
                foc = fun(xoc,varargin{:});
                evals = evals + 1;
                if foc <= fr
                    Vnew = xoc;
                    fnew = foc;
                    how = 'outside contr.';
                else
                    shrink = 1;
                end
            end
            if fr >= f(n+1)
                xic = centroid - paramsNM(3)*(xr - centroid); %inside contraction
                fic = fun(xic,varargin{:});
                evals = evals + 1;
                if fic < f(n+1)
                    Vnew = xic;
                    fnew = fic;
                    how = 'inside contr.';
                else
                    shrink = 1;
                end
            end
        end
    end    
    if shrink == 1
        for j = 2:n+1
            V(j,:) = V(1,:) + paramsNM(4)*(V(j,:) - V(1,:));   %shrink

            f(j) = fun(V(j,:),varargin{:});
        end
        evals = evals + n;            
        how = 'shrink';
    else
        V(n+1,:) = Vnew;
        f(n+1) = fnew;
    end
    
    iterations = iterations+1;
    
    [f, I] = sortrows(f);
    V = V(I,:);

    if display == 2
        fprintf('%6d\t%6d\t%9.3f\t',iterations,evals,f(1));
        for j = 1:n            
            fprintf('%9.3f\t',V(1,j));
        end
        fprintf('%s\n',how);
    end
        
end

x = V(1,:);
fval = f(1);

if iterations < maxIter && evals < maxFunEvals
    exitflag = logical(true);
    output.message = sprintf('Search converged successfully. TolX = %e. TolFun = %e.',tolX,tolFun);
    if display > 1
        fprintf('\n');
        disp(output.message);
    end
else
    output.message = sprintf('Search did not converge: ');
    if iterations >= maxIter
        output.message = [output.message sprintf('maxIter (%d) reached.',maxIter)];
    else
        output.message = [output.message sprintf('maxFunEvals (%d) reached.',maxFunEvals)];
    end    
    if display > 0
        fprintf('\n');
        disp(output.message);
    end
end

if display > 1
    fprintf('\nIterations: %d\nFunction evaluations: %d\nFunction value: %.8g\nParameter values:', iterations, evals, fval);
    for j = 1:n            
        fprintf('\t%.6g',V(1,j));
    end
    fprintf('\n');
end

output.iterations = iterations;
output.funcCount = evals;
output.algorithm = 'Nelder-Mead simplex direct search';

end