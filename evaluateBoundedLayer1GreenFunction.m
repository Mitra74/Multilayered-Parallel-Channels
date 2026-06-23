% evaluateBoundedLayer1GreenFunction.m
%
%   (Renamed from: evaluate3DThreeLayerG.m)
%
% Evaluates the Layer 1 Green's function for the bounded three-layer
% molecular diffusion channel (reflective outer boundary in Layer 1,
% unbounded Layer 3), by numerically inverting the zero-order Hankel
% transform over the supplied lambda grid.
%
% Inputs:
%   lambda_vals  - Hankel transform integration variable
%   A1_vals      - spectral coefficient A1(lambda)
%   B1_vals      - spectral coefficient B1(lambda) (reflective-boundary term)
%   s1_vals      - s1(lambda, omega) for each lambda
%   rho          - radial distance from the transmitter in the (x,y) plane
%   z_eval       - observation depth/height
%   z0           - transmitter depth (z0 <= 0, within Layer 1)
%   alpha1       - Layer 1 diffusion coefficient [m^2/s]
%
% Output:
%   G_val - Green's function value at the requested observation point

function G_val = evaluateBoundedLayer1GreenFunction(lambda_vals, A1_vals, B1_vals, s1_vals, rho, z_eval, z0, alpha1)

n_lambda = length(lambda_vals);
integrand = zeros(n_lambda, 1, 'like', 1i*ones(1));

for i = 1:n_lambda
    lambda = lambda_vals(i);
    A1 = A1_vals(i);
    B1 = B1_vals(i);
    s1 = s1_vals(i);

    % Free-space (fundamental) term centered at the transmitter
    direct_term = exp(-s1*abs(z_eval - z0)) / s1;

    % Boundary-correction terms introduced by the reflective outer
    % boundary and the Layer 1/2 interface
    reflected_term = 4*pi*alpha1*A1*exp(s1*z_eval) + 4*pi*alpha1*B1*exp(-s1*z_eval);

    % Spectral-domain Green's function
    G_lambda = (1/(4*pi*alpha1)) * (direct_term + reflected_term);

    % Hankel transform integrand (zero-order Bessel function)
    integrand(i) = G_lambda * lambda * besselj(0, lambda * rho);
end

% Numerically invert the Hankel transform (trapezoidal rule on a
% uniform lambda grid)
if n_lambda > 1
    dlambda = lambda_vals(2) - lambda_vals(1);
    G_val = dlambda * sum(integrand);
else
    G_val = integrand(1);
end
end
