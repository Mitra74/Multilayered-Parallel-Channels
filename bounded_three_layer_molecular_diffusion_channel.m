% bounded_three_layer_molecular_diffusion_channel.m
%
%   (This file was previously named with a header claiming "unbounded
%    first and third layers, verified with PBS". That description was
%    left over from an earlier version of the script and does not match
%    the solver actually being called here. This script in fact uses
%    the BOUNDED-Layer-1 solver set
%    [solveBoundedThreeLayerGreenFunction.m and its helper functions],
%    and has not yet been validated against particle-based simulation
%    (PBS). The header below has been corrected accordingly.)
%
% Three-layer molecular communication (MC) channel: Green's function /
% impulse response in the time domain.
%
% Configuration: Layer 1 has a reflective (zero-flux) outer boundary;
% Layer 2 has finite thickness L; Layer 3 is semi-infinite (unbounded).
% The point-source transmitter is located in Layer 1.
%
%   Layer 1 (z <= 0)      : bounded source/transmitter layer -> Green's function G
%   Layer 2 (0 < z <= L)  : middle layer                     -> impulse response H2
%   Layer 3 (z > L)       : unbounded outer layer              -> impulse response H3
%
% -------------------------------------------------------------------
% Relation to the oceanic case study in the accompanying paper
% -------------------------------------------------------------------
% The paper's reflective-interface ocean model (Sec. III-B) uses only
% TWO layers: the sea surface microlayer (SML), bounded above by a
% reflective air-water interface at z = 0 and bounded below by the
% interface to the bulk ocean at z = -L_SML; and the bulk ocean itself,
% which is unbounded below that interface. The transmitter sits inside
% the SML, next to the reflective boundary.
%
% This script implements the same underlying physics - one layer with a
% reflective outer boundary next to the transmitter, and one unbounded
% layer far from it - but generalises it to THREE layers and shifts the
% coordinate origin:
%   - Layer 1 here plays the role of the paper's SML: it is bounded by
%     a reflective (zero-flux) boundary at its far edge (z = -L2) and
%     opens onto the next layer at its near edge (z = 0).
%   - Layer 2 is an additional middle layer that is NOT present in the
%     paper's two-layer reflective case; it is included here to
%     demonstrate the general multi-layer formulation.
%   - Layer 3 plays the role of the paper's unbounded bulk ocean.
% The coordinate origin is placed at the open Layer 1/2 interface
% rather than at the reflective boundary, and z increases away from the
% reflective wall instead of decreasing - i.e. the paper's geometry is
% recovered by setting Layer 2 to zero thickness and flipping the sign
% of z. The numerical results (delay, peak concentration, decay
% behaviour, etc.) are equivalent under this relabeling.

clear; clc;

%% Layer diffusion coefficients and degradation rates
D = 10e-9;          % Reference molecular diffusion coefficient [m^2/s]

% Layer 1
Deff_L1 = 1e-9;% 5e-10, 3e-10, 5e-10, 4.61e-9, 5.56e-14;  % Effective diffusion coefficient [m^2/s]
kd_L1 = 0;                                                % Molecular degradation rate [1/s]

% Layer 2
Deff_L2 = 1e-9;% 3e-10, 1e-9, 4.61e-9, 4.39e-13;          % Effective diffusion coefficient [m^2/s]
kd_L2 = 0;                                                % Molecular degradation rate [1/s]

% Layer 3
Deff_L3 = 1e-9;% 2.3e-8, 8.89e-13;                        % Effective diffusion coefficient [m^2/s]
kd_L3 = 0;                                                % Molecular degradation rate [1/s]

%% Channel parameters
params = struct();
params.alpha1 = Deff_L1;   % Layer 1 diffusion coefficient [m^2/s] (z <= 0)
params.alpha2 = Deff_L2;   % Layer 2 diffusion coefficient [m^2/s] (0 < z <= L)
params.alpha3 = Deff_L3;   % Layer 3 diffusion coefficient [m^2/s] (z > L)
params.k1 = params.alpha1; % Layer 1 diffusion coefficient (interface flux terms)
params.k2 = params.alpha2; % Layer 2 diffusion coefficient (interface flux terms)
params.k3 = params.alpha3; % Layer 3 diffusion coefficient (interface flux terms)
% NOTE: k_i is kept separate from alpha_i because this solver was adapted
% from a heat-conduction template (alpha = diffusivity for the PDE,
% k = conductivity for the flux law - two different quantities there).
% For ordinary molecular diffusion there is only one coefficient, so
% k_i = alpha_i is physically required, not a shortcut. This makes
% params.k1/k2/k3 redundant - alpha1/alpha2/alpha3 could be substituted
% directly into the flux-continuity rows with no change to the results.
params.L = 200e-6; %200e-6;  % Layer 1(SML) and 2 thickness [m]

% Interface coupling coefficients (set to 1 for perfect interface coupling)
params.K12 = 1/sqrt(params.alpha1/params.alpha2);  % Interface coupling at z = 0
params.K23 = 1/sqrt(params.alpha2/params.alpha3);  % Interface coupling at z = L

%% Frequency-domain sampling
params.ww = 0.00001;               % Frequency step [rad/s]
params.fw = 10;                    % Maximum frequency [rad/s]
omega_values = params.ww:params.ww:params.fw;
n_frequencies = length(omega_values);

%% Hankel transform discretization
lambda_max = 100000;
n_lambda = 100;
lambda_values = linspace(0.0001, lambda_max, n_lambda);  % Start from small non-zero value

%% Observation points (one per layer)
% Layer 1 (z <= 0) - bounded source/transmitter layer
% Please plot this on paper to see its relation with the configuration in
% the paper. They are identical just the z=0 coordinate is changed.
eval_layer1 = struct();
eval_layer1.x = 0e-6;     eval_layer1.y = 0;        eval_layer1.z =-50e-6; 
eval_layer1.x0 = 0;       eval_layer1.y0 = 0;       eval_layer1.z0 = -190e-6;   % Transmitter location, layer 1

% Layer 2 (0 < z <= L)
eval_layer2 = struct();
eval_layer2.x = 0e-6;     eval_layer2.y = 0;        eval_layer2.z = 50e-6;
eval_layer2.x0 = 0;       eval_layer2.y0 = 0;       eval_layer2.z0 = -190e-6;   % Transmitter location, layer 1

% Layer 3 (z > L) - unbounded outer layer
eval_layer3 = struct();
eval_layer3.x = 0e-6;     eval_layer3.y = 0;        eval_layer3.z = 150e-6;
eval_layer3.x0 = 0;       eval_layer3.y0 = 0;       eval_layer3.z0 = -190e-6;   % Transmitter location, layer 1

%% Frequency-domain Green's function / impulse response
G_omega = zeros(n_frequencies, 1, 'like', 1i*ones(1));    % Layer 1 (bounded source/transmitter layer)
H2_omega = zeros(n_frequencies, 1, 'like', 1i*ones(1));   % Layer 2 (middle layer)
H3_omega = zeros(n_frequencies, 1, 'like', 1i*ones(1));   % Layer 3 (unbounded outer layer)

tic;
for i = 1:n_frequencies
    omega = omega_values(i)*pi;

    try
        results_G = solveBoundedThreeLayerGreenFunction(params, omega, lambda_values, eval_layer1);
        G_omega(i) = results_G.G_val;

        results_H2 = solveBoundedThreeLayerGreenFunction(params, omega, lambda_values, eval_layer2);
        H2_omega(i) = results_H2.H2_val;

        results_H3 = solveBoundedThreeLayerGreenFunction(params, omega, lambda_values, eval_layer3);
        H3_omega(i) = results_H3.H3_val;

    catch ME
        warning('Failed at frequency %.4f rad/s: %s', omega, ME.message);
        G_omega(i) = NaN;
        H2_omega(i) = NaN;
        H3_omega(i) = NaN;
    end
end

total_time = toc;
fprintf('Completed bounded three-layer molecular diffusion channel calculation in %.1f seconds\n', total_time);

%% Time-domain conversion (inverse FFT)
fprintf('Converting to time domain...\n');

domega = omega_values(2) - omega_values(1);  % Frequency step, equal to params.ww

n_points = 3;
eval_z_values = [eval_layer1.z, eval_layer2.z, eval_layer3.z];

G_omega_points = [G_omega.'; H2_omega.'; H3_omega.'];

G_time = zeros(n_points, 2*length(G_omega));  % Pre-allocate for symmetric data

for point_idx = 1:n_points
    % Scale the frequency-domain solution
    C = G_omega_points(point_idx, :) * max(omega_values);

    % Build symmetric frequency-domain data: [0 C conj(fliplr(C))]
    symmetric_data = [0 C conj(fliplr(C))];

    % Inverse FFT to obtain the time-domain response
    time_response = real(ifft(symmetric_data));

    if point_idx == 1
        G_time = zeros(n_points, length(time_response));
    end
    G_time(point_idx, :) = time_response;
end

% Time axis
t_values = 0:1/max(omega_values):1/domega*2+domega;
t_values = t_values(1:size(G_time,2));

fprintf('Time-domain parameters:\n');
fprintf('- domega (frequency step): %.5f rad/s\n', domega);
fprintf('- max(omega): %.1f rad/s\n', max(omega_values));
fprintf('- Time step: %.6f s\n', 1/max(omega_values));
fprintf('- Time points: %d\n', length(t_values));
fprintf('- Time range: %.6f to %.6f s\n', t_values(1), t_values(end));

%% Combined three-layer channel response
t_display = t_values;
G_display = G_time;

layer_names = {'Layer 1 - Green''s function G (bounded source/transmitter layer)', ...
                'Layer 2 - Impulse response H_2 (middle layer)', ...
                'Layer 3 - Impulse response H_3 (unbounded outer layer)'};
layer_colors = {'b-', 'g-', 'r-'};

figure('Name', 'Bounded Three-Layer Molecular Diffusion Channel - Combined Response');
hold on;

for idx = 1:n_points
    z_val = eval_z_values(idx);
    if idx == 1
        label = sprintf('G: z = %.0f \\mum (bounded source/transmitter layer)', z_val*1e6);
    elseif idx == 2
        label = sprintf('H_2: z = +%.0f \\mum (middle layer)', z_val*1e6);
    else
        label = sprintf('H_3: z = +%.0f \\mum (unbounded outer layer)', z_val*1e6);
    end
    plot(t_display, G_display(idx, :), layer_colors{idx}, 'LineWidth', 3, 'DisplayName', label);
end

xlabel('Time [s]', 'FontSize', 14);
ylabel('Green''s function / impulse response [s/m^3]', 'FontSize', 14);
title('Bounded Three-Layer Molecular Diffusion Channel: G, H_2, H_3 vs. Time', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);

y_limits = ylim;
line([0 0], y_limits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');

%% Individual layer responses
figure('Name', 'Layer 1 Green''s Function - Bounded Source/Transmitter Layer');
plot(t_display, G_display(1, :), 'g:', 'LineWidth', 3);
xlabel('Time [s]', 'FontSize', 14);
ylabel('G(r|r_0; t) [s/m^3]', 'FontSize', 14);
title(sprintf('Layer 1 Green''s Function G at z = %.0f \\mum (Bounded Source/Transmitter Layer)', eval_z_values(1)*1e6), 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

figure('Name', 'Layer 2 Impulse Response - Middle Layer');
plot(t_display, G_display(2, :), 'b:', 'LineWidth', 3);
xlabel('Time [s]', 'FontSize', 14);
ylabel('H_2(r; t) [s/m^3]', 'FontSize', 14);
title(sprintf('Layer 2 Impulse Response H_2 at z = +%.0f \\mum (Middle Layer)', eval_z_values(2)*1e6), 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

figure('Name', 'Layer 3 Impulse Response - Unbounded Outer Layer');
plot(t_display, G_display(3, :), 'g:', 'LineWidth', 3);
xlabel('Time [s]', 'FontSize', 14);
ylabel('H_3(r; t) [s/m^3]', 'FontSize', 14);
title(sprintf('Layer 3 Impulse Response H_3 at z = +%.0f \\mum (Unbounded Outer Layer)', eval_z_values(3)*1e6), 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

%% Summary
fprintf('\nBounded Three-Layer Molecular Diffusion Channel - Time-Domain Summary\n');
fprintf('Transmitter location: z0 = %.0f um (Layer 1)\n', eval_layer1.z0*1e6);
fprintf('Layer boundaries: z = 0 and z = %.0f um\n', params.L*1e6);
fprintf('========================================\n');

for point_idx = 1:n_points
    [peak_val, peak_idx] = max(abs(G_time(point_idx, :)));
    peak_time = t_values(peak_idx);
    z_val = eval_z_values(point_idx);

    if point_idx == 1
        fprintf('Layer 1 (G) at z = %.0f um: peak = %.3e at t = %.3f ms\n', z_val*1e6, peak_val, peak_time*1000);
    elseif point_idx == 2
        fprintf('Layer 2 (H2) at z = +%.0f um: peak = %.3e at t = %.3f ms\n', z_val*1e6, peak_val, peak_time*1000);
    else
        fprintf('Layer 3 (H3) at z = +%.0f um: peak = %.3e at t = %.3f ms\n', z_val*1e6, peak_val, peak_time*1000);
    end
end

%% Store results
results = struct();
results.t_values = t_values;
results.G_time = G_time;
results.omega_values = omega_values;
results.G_omega = G_omega;
results.H2_omega = H2_omega;
results.H3_omega = H3_omega;
results.params = params;
results.eval_layer1 = eval_layer1;
results.eval_layer2 = eval_layer2;
results.eval_layer3 = eval_layer3;
results.eval_z_values = eval_z_values;
results.layer_names = layer_names;

fprintf('\nResults stored in the ''results'' variable.\n');
fprintf('Done.\n');
