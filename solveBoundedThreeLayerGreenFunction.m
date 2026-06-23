% solveBoundedThreeLayerGreenFunction.m
%
%   (Renamed from: solve3DThreeLayerGreenFunction.m)
%
% Three-layer molecular diffusion channel solver, bounded configuration.
%
% Configuration: Layer 1 has a reflective (zero-flux) outer boundary at
% z = -L_sml (L_sml=L, Thickness of the SML layer); Layer 2 has finite thickness L; Layer 3 is semi-infinite
% (unbounded). The point-source transmitter is located in Layer 1 -
% for example, the sea surface microlayer (SML) in the ocean case study.
%
% Extensibility:
%   - Other outer-boundary conditions (e.g. a permeable or fixed-flux
%     boundary instead of the reflective one used here) can be added by
%     replacing the boundary-condition row at z = -L_sml with the
%     appropriate equation, and updating evaluateBoundedLayer1GreenFunction.m
%     to match.
%   - The observation points can be placed arbitrarily
%     within their respective layers; only the observation coordinates need to be changed.
%   - For changing the transmitter location, you should make sure that the
%     source containing layer is GREENS FUNCTION, G
%   - Additional layers can be incorporated by extending the linear
%     system below with one field-continuity and one flux-continuity
%     row per extra interface, following the same pattern used for the
%     existing Layer 1/2 and Layer 2/3 interfaces.

function results = solveBoundedThreeLayerGreenFunction(params, omega, lambda_values, evaluation_points)
% Computes the 3D Green's function / impulse response of the bounded
% three-layer channel using the zero-order Hankel transform.
%
% Inputs:
%   params              - struct with the channel parameters:
%                            .alpha1, .alpha2, .alpha3  diffusion coefficients [m^2/s]
%                            .k1, .k2, .k3              diffusion coefficients used in
%                                                        interface flux-matching terms
%                            .L                          Layer 2 thickness [m]
%                            .K12, .K23                 interface coupling coefficients
%   omega                - angular frequency [rad/s]
%   lambda_values        - Hankel transform integration grid
%   evaluation_points     - struct with observation coordinates (x,y,z)
%                           and transmitter coordinates (x0,y0,z0)
%
% Output:
%   results - struct with the spectral coefficients (A1, B1, A2, B2, B3),
%             the evaluated Green's function / impulse response at the
%             requested observation point, and bookkeeping fields
%             (region, rho, params, omega).

%% Layer parameters
alpha1 = params.alpha1;  % Layer 1 diffusion coefficient [m^2/s],-L_sml< z <= 0
alpha2 = params.alpha2;  % Layer 2 diffusion coefficient [m^2/s], 0 < z <= L
alpha3 = params.alpha3;  % Layer 3 diffusion coefficient [m^2/s], z > L
k1 = params.k1;          % Layer 1 diffusion coefficient (interface flux terms)
k2 = params.k2;          % Layer 2 diffusion coefficient (interface flux terms)
k3 = params.k3;          % Layer 3 diffusion coefficient (interface flux terms)
L = params.L;            % Layer 2 thickness [m]
K12 = params.K12;        % Interface coupling coefficient at z = 0
K23 = params.K23;        % Interface coupling coefficient at z = L

% Observation and transmitter coordinates
x_eval = evaluation_points.x;
y_eval = evaluation_points.y;
z_eval = evaluation_points.z;
x0 = evaluation_points.x0;
y0 = evaluation_points.y0;
z0 = evaluation_points.z0;

% Radial distance between the observation point and the transmitter
rho = sqrt((x_eval - x0)^2 + (y_eval - y0)^2);

if z0 > 0
    error('Transmitter must be located in layer 1: z0 <= 0 (z0 = %.3f)', z0);
end

%% Spectral coefficients (one set per Hankel transform variable lambda)
n_lambda = length(lambda_values);
A1_lambda = zeros(n_lambda, 1, 'like', 1i*ones(1));
B1_lambda = zeros(n_lambda, 1, 'like', 1i*ones(1));  % Reflective-boundary coefficient, Layer 1
A2_lambda = zeros(n_lambda, 1, 'like', 1i*ones(1));
B2_lambda = zeros(n_lambda, 1, 'like', 1i*ones(1));
B3_lambda = zeros(n_lambda, 1, 'like', 1i*ones(1));

% s_p(lambda, omega) = sqrt(lambda^2 + i*omega/alpha_p) for each layer p
s1_vals = sqrt(lambda_values.^2 + 1i*omega/alpha1);
s2_vals = sqrt(lambda_values.^2 + 1i*omega/alpha2);
s3_vals = sqrt(lambda_values.^2 + 1i*omega/alpha3);

% Outer-boundary location for the reflective condition in Layer 1
% (equal to params.L for the configuration used in this script)
L_sml = 200e-6; 

%% Solve the 5x5 boundary-condition system for each lambda
for i = 1:n_lambda
    lambda = lambda_values(i);
    s1 = s1_vals(i);
    s2 = s2_vals(i);
    s3 = s3_vals(i);

    A_matrix = zeros(5, 5, 'like', 1i*ones(1));
    b_vector = zeros(5, 1, 'like', 1i*ones(1));

    % Row 1 - reflective (zero-flux) boundary at z = -L_sml: dG1/dz = 0
    A_matrix(1, 1) = 4*pi*alpha1*s1*exp(-s1*L_sml);
    A_matrix(1, 2) = -4*pi*alpha1*s1*exp(s1*L_sml);
    A_matrix(1, 3) = 0;
    A_matrix(1, 4) = 0;
    A_matrix(1, 5) = 0;
    b_vector(1) = -exp(-s1*abs(-z0-L_sml));

    % Row 2 - concentration continuity at z = 0: G1(0-) = K12 * H2(0+)
    A_matrix(2, 1) = 4*pi*alpha1;
    A_matrix(2, 2) = 4*pi*alpha1;
    A_matrix(2, 3) = -K12*4*pi*alpha1;
    A_matrix(2, 4) = -K12*4*pi*alpha1;
    A_matrix(2, 5) = 0;
    b_vector(2) = -exp(-s1*abs(z0))/s1;

    % Row 3 - flux continuity at z = 0: k1*dG1/dz(0-) = k2*dH2/dz(0+)
    A_matrix(3, 1) = 4*pi*alpha1*k1*s1;
    A_matrix(3, 2) = -4*pi*alpha1*k1*s1;
    A_matrix(3, 3) = k2*4*pi*alpha1*s2;
    A_matrix(3, 4) = -k2*4*pi*alpha1*s2;
    A_matrix(3, 5) = 0;
    b_vector(3) = k1*exp(-s1*abs(z0));

    % Row 4 - concentration continuity at z = L: H2(L-) = K23 * H3(L+)
    A_matrix(4, 1) = 0;
    A_matrix(4, 2) = 0;
    A_matrix(4, 3) = exp(-s2*L);
    A_matrix(4, 4) = exp(s2*L);
    A_matrix(4, 5) = -K23;
    b_vector(4) = 0;

    % Row 5 - flux continuity at z = L: k2*dH2/dz(L-) = k3*dH3/dz(L+)
    A_matrix(5, 1) = 0;
    A_matrix(5, 2) = 0;
    A_matrix(5, 3) = -s2*k2*exp(-s2*L);
    A_matrix(5, 4) = s2*k2*exp(s2*L);
    A_matrix(5, 5) = s3*k3;
    b_vector(5) = 0;

    try
        coef = pinv(A_matrix) * b_vector;
        A1_lambda(i) = coef(1);
        B1_lambda(i) = coef(2);
        A2_lambda(i) = coef(3);
        B2_lambda(i) = coef(4);
        B3_lambda(i) = coef(5);
    catch ME
        warning('Failed to solve for lambda = %.4f: %s', lambda, ME.message);
        A1_lambda(i) = NaN;
        B1_lambda(i) = NaN;
        A2_lambda(i) = NaN;
        B2_lambda(i) = NaN;
        B3_lambda(i) = NaN;
    end
end

%% Evaluate the Green's function / impulse response in the requested layer
valid_idx = ~isnan(A1_lambda) & ~isnan(B1_lambda) & ~isnan(A2_lambda) & ~isnan(B2_lambda) & ~isnan(B3_lambda);
lambda_valid = lambda_values(valid_idx);

if z_eval <= 0
    G_val = evaluateBoundedLayer1GreenFunction(lambda_valid, A1_lambda(valid_idx), B1_lambda(valid_idx), ...
                                               s1_vals(valid_idx), rho, z_eval, z0, alpha1);
    H2_val = NaN;
    H3_val = NaN;
    region = 'Layer 1 (z <= 0) - source/transmitter layer';
elseif z_eval > 0 && z_eval <= L
    H2_val = evaluateBoundedLayer2ImpulseResponse(lambda_valid, A2_lambda(valid_idx), B2_lambda(valid_idx), ...
                                                  s2_vals(valid_idx), rho, z_eval);
    G_val = NaN;
    H3_val = NaN;
    region = 'Layer 2 (0 < z <= L)';
else
    H3_val = evaluateBoundedLayer3ImpulseResponse(lambda_valid, B3_lambda(valid_idx), ...
                                                  s3_vals(valid_idx), rho, z_eval, L);
    G_val = NaN;
    H2_val = NaN;
    region = 'Layer 3 (z > L)';
end

%% Package results
results = struct();
results.lambda_values = lambda_values;
results.A1_lambda = A1_lambda;
results.B1_lambda = B1_lambda;
results.A2_lambda = A2_lambda;
results.B2_lambda = B2_lambda;
results.B3_lambda = B3_lambda;
results.G_val = G_val;
results.H2_val = H2_val;
results.H3_val = H3_val;
results.region = region;
results.rho = rho;
results.params = params;
results.omega = omega;
end
