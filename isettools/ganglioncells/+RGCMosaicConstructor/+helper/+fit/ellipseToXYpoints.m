%
% RGCMosaicConstrucor.helper.fit.ellipseToXYpoints(xyPoints, varargin)
% 
function [z, a, b, alpha] = ellipseToXYpoints(xyPoints, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('maxIterations', 200, @isnumeric);
    p.addParameter('tolerance', 1e-5, @isnumeric);
    p.addParameter('constraint', 'bookstein', @(x)(ismember(x, {'bookstein', 'trace'})));
    p.addParameter('nonLinear', true, @islogical);
    p.parse(varargin{:});

	rows = size(xyPoints,1);
	cols = size(xyPoints,2);
	assert(rows == 2, 'xyPoints must be a 2 x N matrix');
	assert(cols >=6, 'xyPoints must have at least 6 points');

	fitParams = struct(...
		'nonLinear', p.Results.nonLinear, ...
		'constraint', p.Results.constraint, ...
		'maxIterations', p.Results.maxIterations, ...
		'tolerance', p.Results.tolerance);

	% Remove centroid
	centroid = mean(xyPoints, 2);
	xyPoints = bsxfun(@minus, xyPoints, centroid);

	% Obtain a linear estimate
	switch fitParams.constraint
    	case 'bookstein'
        	[z, a, b, alpha] = fitbookstein(xyPoints);
    	case 'trace'
       		[z, a, b, alpha] = fitggk(xyPoints);
	end % switch


	if (fitParams.nonLinear)
		% Initial conditions
	    z0     = z;
	    a0     = a;
	    b0     = b;
	    alpha0 = alpha;

	    % Fit
    	[z, a, b, alpha, converged, isCircle] = fitNonLinear(xyPoints, z0, a0, b0, alpha0, params);

    	% Return linear estimate if GN doesn't converge or if the data points fall on a circle
	    if (~converged) || (isCircle)
	        fprintf('*** FailureToConverge: Gauss-Newton did not converge, returning linear estimate.\n');
	        z = z0;
	        a = a0;
	        b = b0;
	        alpha = alpha0;
	    end
	end % if (fitParams.nonLinear)

	% Add the centroid back on
	z = z + centroid;
end


function [z, a, b, alpha] = fitbookstein(x)
	%FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
	%   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A

	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Define the coefficient matrix B, such that we solve the system
	% B *[v; w] = 0, with the constraint norm(w) == 1
	B = [x1, x2, ones(m, 1), x1.^2, sqrt(2) * x1 .* x2, x2.^2];

	% To enforce the constraint, we need to take the QR decomposition
	[Q, R] = qr(B);

	% Decompose R into blocks
	R11 = R(1:3, 1:3);
	R12 = R(1:3, 4:6);
	R22 = R(4:6, 4:6);

	% Solve R22 * w = 0 subject to norm(w) == 1
	[U, S, V] = svd(R22);
	w = V(:, 3);

	% Solve for the remaining variables
	v = -R11 \ R12 * w;

	% Fill in the quadratic form
	A        = zeros(2);
	A(1)     = w(1);
	A([2 3]) = 1 / sqrt(2) * w(2);
	A(4)     = w(3);
	bv       = v(1:2);
	c        = v(3);

	% Find the parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end

function [z, a, b, alpha] = fitggk(x)
	% Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Coefficient matrix
	B = [2 * x1 .* x2, x2.^2 - x1.^2, x1, x2, ones(m, 1)];
	v = B \ -x1.^2;

	% For clarity, fill in the quadratic form variables
	A        = zeros(2);
	A(1,1)   = 1 - v(2);
	A([2 3]) = v(1);
	A(2,2)   = v(2);
	bv       = v(3:4);
	c        = v(5);

	% find parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end


function [z, a, b, alpha, converged, isCircle] = fitNonLinear(x, z0, a0, b0, alpha0, params)
	% Gauss-Newton least squares ellipse fit minimising geometric distance 

	% Get initial rotation matrix
	Q0 = [cos(alpha0), -sin(alpha0); sin(alpha0) cos(alpha0)];
	m = size(x, 2);

	% Get initial phase estimates
	phi0 = angle( [1 i] * Q0' * (x - repmat(z0, 1, m)) )';
	u = [phi0; alpha0; a0; b0; z0];

	% Iterate using Gauss Newton
	converged = false;

	for nIts = 1:params.maxIterations
	    % Find the function and Jacobian
	    [f, J, isCircle] = computeJacobian(u);
    
    	if (isCircle)
    		fprintf('Ellipse is near-circular - nonlinear fit may not succeed\n.')
    	end

	    % Solve for the step and update u
	    h = -J \ f;
	    u = u + h;
    
	    % Check for convergence
	    delta = norm(h, inf) / norm(u, inf);
	    if delta < params.tolerance
	        converged = true;
	        break
	    end
	end

	alpha = u(end-4);
	a = u(end-3);
	b = u(end-2);
	z = u(end-1:end);

	% ---- Nested function ---
	function [f, J, isCircle] = computeJacobian(u)
        % Define the system of nonlinear equations and Jacobian. 

        % Tolerance for whether it is a circle
        circTol = 1e-5;
        
        % Unpack parameters from u
        phi   = u(1:end-5);
        alpha = u(end-4);
        a     = u(end-3);
        b     = u(end-2);
        z     = u(end-1:end);
        
        % If it is a circle, the Jacobian will be singular, and the
        % Gauss-Newton step won't work. 
        %TODO: This can be fixed by switching to a Levenberg-Marquardt
        %solver
        if (abs(a - b) / (a + b) < circTol)
            isCircle = true;
        else
        	isCircle = false;
        end

        % Convenience trig variables
        c = cos(phi);
        s = sin(phi);
        ca = cos(alpha);
        sa = sin(alpha);
        
        % Rotation matrices
        Q    = [ca, -sa; sa, ca];
        Qdot = [-sa, -ca; ca, -sa];

        % Preallocate function and Jacobian variables
        f = zeros(2 * m, 1);
        J = zeros(2 * m, m + 5);
        for i = 1:m
            rows = (2*i-1):(2*i);
            % Equation system - vector difference between point on ellipse
            % and data point
            f((2*i-1):(2*i)) = x(:, i) - z - Q * [a * cos(phi(i)); b * sin(phi(i))];
            
            % Jacobian
            J(rows, i) = -Q * [-a * s(i); b * c(i)];
            J(rows, (end-4:end)) = ...
                [-Qdot*[a*c(i); b*s(i)], -Q*[c(i); 0], -Q*[0; s(i)], [-1 0; 0 -1]];
        end
    end % ---- Nested function ---
end

function [z, a, b, alpha] = conic2parametric(A, bv, c)
	% Diagonalise A - find Q, D such at A = Q' * D * Q
	[Q, D] = eig(A);
	Q = Q';

	% If the determinant < 0, it's not an ellipse
	if prod(diag(D)) <= 0 
	    error('NotEllipse', 'Linear fit did not produce an ellipse');
	end

	% We have b_h' = 2 * t' * A + b'
	t = -0.5 * (A \ bv);

	c_h = t' * A * t + bv' * t + c;

	z = t;
	a = sqrt(-c_h / D(1,1));
	b = sqrt(-c_h / D(2,2));
	alpha = atan2(Q(1,2), Q(1,1));
end % conic2parametric