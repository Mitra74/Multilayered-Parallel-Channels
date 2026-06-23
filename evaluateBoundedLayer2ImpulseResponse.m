% evaluateBoundedLayer2ImpulseResponse.m
%
%   (Renamed from: evaluate3DThreeLayerH2.m)
%
% Evaluates the Layer 2 impulse response (middle layer, 0 < z <= L) of
% the bounded three-layer molecular diffusion channel, by numerically
% inverting the zero-order Hankel transform over the supplied lambda
% grid.
%
% Inputs:
%   lambda_vals  - Hankel transform integration variable
%   A2_vals      - spectral coefficient A2(lambda)
%   B2_vals      - spectral coefficient B2(lambda)
%   s2_vals      - s2(lambda, omega) for each lambda
%   rho          - radial distance from the transmitter in the (x,y) plane
%   z_eval       - observation depth/height within Layer 2
%
% Output:
%   H2_val - impulse response value at the requested observation point

function H2_val = evaluateBoundedLayer2ImpulseResponse(lambda_vals, A2_vals, B2_vals, s2_vals, rho, z_eval)

n_lambda = length(lambda_vals);
integrand = zeros(n_lambda, 1, 'like', 1i*ones(1));

for i = 1:n_lambda
    lambda = lambda_vals(i);
    A2_val = A2_vals(i);
    B2_val = B2_vals(i);
    s2 = s2_vals(i);

    % Two-term exponential solution within the middle layer
    exp_term = A2_val * exp(-s2 * z_eval) + B2_val * exp(s2 * z_eval);

    % Hankel transform integrand (zero-order Bessel function)
    J0_term = besselj(0, lambda * rho);

    integrand(i) = exp_term * J0_term * lambda;
end

% Numerically invert the Hankel transform
H2_val = trapz(lambda_vals, integrand);
end
