% evaluateBoundedLayer3ImpulseResponse.m
%
%   (Renamed from: evaluate3DThreeLayerH3.m)
%
% Evaluates the Layer 3 impulse response (outer layer, z > L) of the
% bounded three-layer molecular diffusion channel, by numerically
% inverting the zero-order Hankel transform over the supplied lambda
% grid. The exponential decay is referenced to z = L (rather than
% z = 0) so the solution decays as z -> infinity while remaining
% continuous with Layer 2 at the interface z = L.
%
% Inputs:
%   lambda_vals  - Hankel transform integration variable
%   B3_vals      - spectral coefficient B3(lambda)
%   s3_vals      - s3(lambda, omega) for each lambda
%   rho          - radial distance from the transmitter in the (x,y) plane
%   z_eval       - observation depth/height within Layer 3 (z_eval >= L)
%   L            - location of the Layer 2 / Layer 3 interface
%
% Output:
%   H3_val - impulse response value at the requested observation point

function H3_val = evaluateBoundedLayer3ImpulseResponse(lambda_vals, B3_vals, s3_vals, rho, z_eval, L)

n_lambda = length(lambda_vals);
integrand = zeros(n_lambda, 1, 'like', 1i*ones(1));

for i = 1:n_lambda
    lambda = lambda_vals(i);
    B3_val = B3_vals(i);
    s3 = s3_vals(i);

    % Decaying exponential term referenced to the Layer 2/3 interface
    exp_term = B3_val * exp(-s3 * (z_eval - L));

    % Hankel transform integrand (zero-order Bessel function)
    J0_term = besselj(0, lambda * rho);

    integrand(i) = exp_term * J0_term * lambda;
end

% Numerically invert the Hankel transform
H3_val = trapz(lambda_vals, integrand);
end
