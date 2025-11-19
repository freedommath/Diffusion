% Temperature and Interface Comparison: 330C vs 500C (Sharp vs Gradient)
% This script compares 4 cases:
% 1. 330°C Sharp Interface
% 2. 330°C Gradient Interface
% 3. 500°C Sharp Interface
% 4. 500°C Gradient Interface

clear;
clc;
close all;

%% Temperature and Material Parameters
% Temperature settings
T1 = 330 + 273.15;  % 330°C to K (603.15 K)
T2 = 500 + 273.15;  % 500°C to K (773.15 K)
R = 8.314;          % Gas constant (J/mol·K)

% Diffusion parameters for Ni
D0_Ni = 4.7e-9;     % Pre-exponential factor for Ni (m²/s)
Q_Ni = 143000;      % Activation energy for Ni (J/mol)

% Diffusion parameters for Ti
D0_Ti = 2.7e-7;     % Pre-exponential factor for Ti (m²/s)
Q_Ti = 205000;      % Activation energy for Ti (J/mol)

% Diffusion coefficients calculation
DNi_330 = D0_Ni * exp(-Q_Ni / (R * T1));
DTi_330 = D0_Ti * exp(-Q_Ti / (R * T1));
DNi_500 = D0_Ni * exp(-Q_Ni / (R * T2));
DTi_500 = D0_Ti * exp(-Q_Ti / (R * T2));

fprintf('=== Diffusion Coefficient Comparison ===\n');
fprintf('330°C: D_Ni = %.2e m²/s, D_Ti = %.2e m²/s, Ratio = %.1f\n', DNi_330, DTi_330, DNi_330/DTi_330);
fprintf('500°C: D_Ni = %.2e m²/s, D_Ti = %.2e m²/s, Ratio = %.1f\n', DNi_500, DTi_500, DNi_500/DTi_500);

%% Simulation Parameters
total_length = 410e-9;     % Total specimen length (m)
gradient_thickness = 10e-9;% Gradient interface thickness (m)
sim_time = 100;            % Total simulation time (s)

nx = 1011;                 % Number of spatial grid points
dx = total_length / (nx - 1); % Spatial grid spacing (m)
x = linspace(0, total_length, nx); % Spatial grid creation

%% Initial Concentration Profiles
% Sharp Interface initial condition
conc_Ni_sharp_initial = ones(1, nx);
conc_Ni_sharp_initial(x > 200e-9) = 0;

% Gradient Interface initial condition
conc_Ni_gradient_initial = ones(1, nx);
gradient_start_pos = 200e-9;
gradient_end_pos = gradient_start_pos + gradient_thickness;
gradient_indices = x > gradient_start_pos & x <= gradient_end_pos;
conc_Ni_gradient_initial(gradient_indices) = 1 - (x(gradient_indices) - gradient_start_pos) / gradient_thickness;
conc_Ni_gradient_initial(x > gradient_end_pos) = 0;

%% Run Simulations using the created function
% --- 330°C Simulations ---
conc_Ni_330_sharp = run_diffusion_sim(conc_Ni_sharp_initial, DNi_330, DTi_330, sim_time, dx);
conc_Ni_330_gradient = run_diffusion_sim(conc_Ni_gradient_initial, DNi_330, DTi_330, sim_time, dx);

% --- 500°C Simulations ---
conc_Ni_500_sharp = run_diffusion_sim(conc_Ni_sharp_initial, DNi_500, DTi_500, sim_time, dx);
conc_Ni_500_gradient = run_diffusion_sim(conc_Ni_gradient_initial, DNi_500, DTi_500, sim_time, dx);


%% Diffusion Amount Calculation (converted to nm)
diffused_330_sharp = sum(abs(conc_Ni_330_sharp - conc_Ni_sharp_initial)) * dx * 1e9; % nm
diffused_330_gradient = sum(abs(conc_Ni_330_gradient - conc_Ni_gradient_initial)) * dx * 1e9; % nm
diffused_500_sharp = sum(abs(conc_Ni_500_sharp - conc_Ni_sharp_initial)) * dx * 1e9; % nm
diffused_500_gradient = sum(abs(conc_Ni_500_gradient - conc_Ni_gradient_initial)) * dx * 1e9; % nm

%% Results Output (Identical to original)
fprintf('\n--- 330°C Results ---\n');
fprintf('Sharp Interface   : %.2f nm\n', diffused_330_sharp);
fprintf('Gradient Interface: %.2f nm\n', diffused_330_gradient);
fprintf('Sharp/Gradient Ratio: %.2f\n', diffused_330_sharp/diffused_330_gradient);

fprintf('\n--- 500°C Results ---\n');
fprintf('Sharp Interface   : %.2f nm\n', diffused_500_sharp);
fprintf('Gradient Interface: %.2f nm\n', diffused_500_gradient);
fprintf('Sharp/Gradient Ratio: %.2f\n', diffused_500_sharp/diffused_500_gradient);

fprintf('\n--- Temperature Comparison ---\n');
fprintf('Sharp Interface (500°C/330°C):    %.2f\n', diffused_500_sharp/diffused_330_sharp);
fprintf('Gradient Interface (500°C/330°C): %.2f\n', diffused_500_gradient/diffused_330_gradient);

%% Visualization (Identical to original)
figure('Name', 'Temperature Comparison: 330°C vs 500°C Interface Effects', 'Position', [100, 100, 1200, 600]);

% 330°C comparison
subplot(1, 2, 1);
plot(x * 1e9, conc_Ni_sharp_initial, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Initial (2D)');
hold on;
plot(x * 1e9, conc_Ni_gradient_initial, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial (3D)');
plot(x * 1e9, conc_Ni_330_sharp, 'r-', 'LineWidth', 2, 'DisplayName', 'Sharp after 100s');
plot(x * 1e9, conc_Ni_330_gradient, 'g-', 'LineWidth', 2, 'DisplayName', 'Gradient after 100s');
hold off;
title('330°C: 2D vs 3D Interface', 'FontSize', 14);
xlabel('Position (nm)', 'FontSize', 12);
ylabel('Concentration of Ni', 'FontSize', 12);
legend('Location', 'best');
grid on;
xlim([170, 240]);
ylim([-0.1, 1.1]);

% 500°C comparison
subplot(1, 2, 2);
plot(x * 1e9, conc_Ni_sharp_initial, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Initial (2D)');
hold on;
plot(x * 1e9, conc_Ni_gradient_initial, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Initial (3D)');
plot(x * 1e9, conc_Ni_500_sharp, 'r-', 'LineWidth', 2, 'DisplayName', 'Sharp after 100s');
plot(x * 1e9, conc_Ni_500_gradient, 'g-', 'LineWidth', 2, 'DisplayName', 'Gradient after 100s');
hold off;
title('500°C: 2D vs 3D Interface', 'FontSize', 14);
xlabel('Position (nm)', 'FontSize', 12);
ylabel('Concentration of Ni', 'FontSize', 12);
legend('Location', 'best');
grid on;
xlim([170, 240]);
ylim([-0.1, 1.1]);

% Annotation for 330°C
annotation('textbox', [0.33, 0.5, 0.2, 0.1], 'String', ...
    sprintf('2D Interface:\n%.2f nm\n\n3D Interface:\n%.2f nm', ...
    diffused_330_sharp, diffused_330_gradient), ...
    'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'k', ...
    'BackgroundColor', 'w', 'FitBoxToText', 'on');

% Annotation for 500°C
annotation('textbox', [0.77, 0.5, 0.2, 0.1], 'String', ...
    sprintf('2D Interface:\n%.2f nm\n\n3D Interface:\n%.2f nm', ...
    diffused_500_sharp, diffused_500_gradient), ...
    'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'k', ...
    'BackgroundColor', 'w', 'FitBoxToText', 'on');

%% ===== Reusable Simulation Function =====
function conc_final = run_diffusion_sim(conc_initial, D1, D2, time_total, dx)
    % This function runs the 1D diffusion simulation
    % Inputs:
    %   conc_initial: Initial concentration profile (vector)
    %   D1, D2: Diffusion coefficients for species 1 and 2 (e.g., D_Ni, D_Ti)
    %   time_total: Total simulation time (s)
    %   dx: Spatial grid spacing (m)
    % Output:
    %   conc_final: Concentration profile after simulation

    % --- Time Step Calculation for Stability ---
    D_max = max(D1, D2);
    dt = dx^2 / (2 * D_max) * 0.3; % Safety factor of 0.3
    nt = round(time_total / dt);

    conc_current = conc_initial;

    % --- Time Evolution Loop ---
    for t = 1:nt
        % Effective diffusion coefficient (depends on local concentration)
        D_eff = D1 * conc_current + D2 * (1 - conc_current);
        
        % Calculate second derivative using central difference
        d2C_dx2 = (circshift(conc_current, -1) - 2 * conc_current + circshift(conc_current, 1)) / dx^2;
        
        % Apply Neumann (zero-flux) boundary conditions
        d2C_dx2(1) = 2 * (conc_current(2) - conc_current(1)) / dx^2;
        d2C_dx2(end) = 2 * (conc_current(end-1) - conc_current(end)) / dx^2;
        
        % Update concentration using Forward Euler method
        conc_current = conc_current + dt * D_eff .* d2C_dx2;
    end
    
    conc_final = conc_current;
end