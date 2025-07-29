clear all;
clc;

% --- Model Parameters ---
m = 2;          % Mass (kg)
l = 0.2;        % Arm length (m)
g = 9.81;       % Gravity (m/s²)
Ix = 1.2;       % Moment of inertia (N·m²/rad)
Iy = 1.2;       % Moment of inertia (N·m²/rad)
Iz = 2.2;       % Moment of inertia (N·m²/rad)
K1 = 0.01;      % Drag coefficient (N·s/m)
K2 = 0.01;      % Drag coefficient (N·s/m)
K3 = 0.01;      % Drag coefficient (N·s/m)
K4 = 0.012;     % Drag coefficient (N·m·s/rad)
K5 = 0.012;     % Drag coefficient (N·m·s/rad)
K6 = 0.012;     % Drag coefficient (N·m·s/rad)
Jr = 0.2;       % Rotor inertia (N·m²/rad)
Omega_r = 40;    % Total rotor speed (rad/s, simplified)
beta_min_z= sqrt(0.5);
beta_hat_z = sqrt(0.5);

% --- Simulation Parameters ---
Ts = 0.01;      % Sampling time (s)
tspan = 0:Ts:30; % Time span for simulation (0 to 30 seconds)
N = length(tspan); % Number of time steps

% --- Initial State ---
% State vector: X = [x, x_dot, y, y_dot, z, z_dot, phi, phi_dot, theta, theta_dot, psi, psi_dot]'
X0 = zeros(12, 1);

% --- Reference Trajectory ---
% Define a desired trajectory for the quadrotor to follow with variation in both x and y
x_ref = 0.2*sin(0.2*tspan)'; % Modified to have sinusoidal movement in x
y_ref = 0.3*sin(0.15*tspan)'; % Modified to have sinusoidal movement in y
z_ref = ones(N, 1);  % Hover at 1 meter
psi_ref = zeros(N, 1);
phi_ref = zeros(N, 1);
theta_ref = zeros(N, 1);

% Add a step change to demonstrate robustness
change_idx = round(N/2);
z_ref(change_idx:end) = 2;  % Step change to 2 meters at half-time

% --- SMC Parameters ---
% Modified parameters for better response
cz = 1;         % z-direction parameter (increased)
c_psi = 1;      % yaw parameter (increased)
epsilon_1 = 0.8;  % Increased from 0.8
epsilon_2 = 0.8;  % Increased from 0.8
epsilon_3 = 0.5;  % Increased from 0.5
epsilon_4 = 0.5;  % Increased from 0.5
eta_1 = 2;      % Increased from 2
eta_2 = 2;      % Increased from 2
eta_3 = 5;        % Increased from 5
eta_4 = 5;        % Increased from 5ㄈ
c3 = 1;         % constant for phi SMC (increased)
c4 = 6;           % constant for phi SMC (increased)
c7 = 1;         % constant for theta SMC (increased)
c8 = 6;           % constant for theta SMC (increased)

% --- Define the continuous-time nonlinear state-space model ---
f_continuous = @(X, U) [
    X(2);   % dx/dt = x_dot
    (cos(X(7)) * sin(X(9)) * cos(X(11)) + sin(X(7)) * sin(X(11))) * U(1)/m -K1*X(2)/m;  % dx_dot/dt
    X(4);   % dy/dt = y_dot
    (cos(X(7)) * sin(X(9)) * sin(X(11)) - sin(X(7)) * cos(X(11))) * U(1) / m - K2*X(4)/m;  % dy_dot/dt
    X(6);   % dz/dt = z_dot
    (cos(X(7)) * cos(X(9))) * U(1) / m - g - K3 * X(6) / m;  % dz_dot/dt
    X(8);   % dphi/dt = phi_dot
    X(10) * X(12) * (Iy - Iz) / Ix + Jr / Ix * X(10) * Omega_r + l / Ix * U(2) - K4 / Ix * X(8);  % dphi_dot/dt
    X(10);  % dtheta/dt = theta_dot
    X(8) * X(12) * (Iz - Ix) / Iy - Jr / Iy * X(8) * Omega_r + l / Iy * U(3) - K5 / Iy * X(10);  % dtheta_dot/dt
    X(12);  % dpsi/dt = psi_dot
    X(8) * X(10) * (Ix - Iy) / Iz + 1 / Iz * U(4) - K6 / Iz * X(12)  % dpsi_dot/dt
];

% --- Construct Continuous-Time Linear Part (Ac, Bc) ---
Ac = zeros(12, 12);
% Position and velocity relationships
Ac(1, 2) = 1;  % dx/dt = x_dot
Ac(2, 2) = -K1/m;  % Drag term for x_dot
Ac(3, 4) = 1;  % dy/dt = y_dot
Ac(4, 4) = -K2/m;  % Drag term for y_dot
Ac(5, 6) = 1;  % dz/dt = z_dot
Ac(6, 6) = -K3/m;  % Drag term for z_dot
% Attitude and angular velocity relationships
Ac(7, 8) = 1;  % dphi/dt = phi_dot
Ac(8, 8) = -K4/Ix;  % Drag term for phi_dot
Ac(9, 10) = 1;  % dtheta/dt = theta_dot
Ac(10, 10) = -K5/Iy;  % Drag term for theta_dot
Ac(11, 12) = 1;  % dpsi/dt = psi_dot
Ac(12, 12) = -K6/Iz;  % Drag term for psi_dot

Bc = zeros(12, 4);     % u1 (thrust) affects z acceleration (modified - direct effect)
Bc(8, 2) = l/Ix;    % u2 (roll torque) affects phi acceleration
Bc(10, 3) = l/Iy;   % u3 (pitch torque) affects theta acceleration
Bc(12, 4) = 1/Iz;   % u4 (yaw torque) affects psi acceleration

% --- Discretize using c2d (First-Order Hold method) ---
sys_continuous = ss(Ac, Bc, eye(12), zeros(12, 4));
sys_discrete = c2d(sys_continuous, Ts, 'zoh');
Ad = sys_discrete.A;
Bd = sys_discrete.B;

% --- Define the nonlinear part ---
f_nonlinear = @(X, U) [
    0;
    (cos(X(7)) * sin(X(9)) * cos(X(11)) + sin(X(7)) * sin(X(11))) * U(1) / m;
    0;
    (cos(X(7)) * sin(X(9)) * sin(X(11)) - sin(X(7)) * cos(X(11))) * U(1) / m;
    0;
    (cos(X(7)) * cos(X(9))) * U(1) / m - g;
    0;
    X(10) * X(12) * (Iy - Iz) / Ix + Jr / Ix * X(10) * Omega_r;
    0;
    X(8) * X(12) * (Iz - Ix) / Iy - Jr / Iy * X(8) * Omega_r;
    0;
    X(8) * X(10) * (Ix - Iy) / Iz
];

% --- Discrete-time model using c2d ---
discreteC2D = @(X, U, dt) Ad * X + Bd * U + dt * f_nonlinear(X, U);

% --- Simulation with SMC ---
X_history = zeros(12, N);
U_history = zeros(4, N-1);  % Adjusted size - control inputs are for N-1 steps
X_history(:, 1) = X0;
X = X0;

% Arrays to record history of desired angles and virtual control inputs
theta_des_history = zeros(1, N-1);
phi_des_history = zeros(1, N-1);
virtual_control_history = zeros(4, N-1); % [x_dot_dot_des, y_dot_dot_des, z_dot_dot_des, psi_dot_dot_des]

% Create array to store sliding mode coefficients
c_history = zeros(N-1, 8); % Store c1 to c8 history

% --- Simulation loop with updated SMC parameters ---
for k = 1:N-1
    % Get current state
    x_k = X(1); x_dot_k = X(2);
    y_k = X(3); y_dot_k = X(4);
    z_k = X(5); z_dot_k = X(6);
    phi_k = X(7); phi_dot_k = X(8);
    theta_k = X(9); theta_dot_k = X(10);
    psi_k = X(11); psi_dot_k = X(12);
    
    % Reference trajectory and reference derivatives calculation
    x_ref_k = x_ref(k);
    y_ref_k = y_ref(k);
    z_ref_k = z_ref(k);
    psi_ref_k = psi_ref(k);
    
    % Calculate reference velocities (now non-zero due to sinusoidal references)
    if k > 1
        x_ref_dot_k = (x_ref(k) - x_ref(k-1)) / Ts;
        y_ref_dot_k = (y_ref(k) - y_ref(k-1)) / Ts;
    else
        x_ref_dot_k = 0;
        y_ref_dot_k = 0;
    end
    z_ref_dot_k = 0; % Still assuming zero velocity for z reference
    psi_ref_dot_k = 0;
    theta_ref_dot_k = 0;
    phi_ref_dot_k = 0;
    % Calculate tracking errors
    e_x = x_k - x_ref_k;
    e_y = y_k - y_ref_k;
    e_z = z_k - z_ref_k;
    e_psi = psi_k - psi_ref_k;
    
    e_x_dot = x_dot_k - x_ref_dot_k;
    e_y_dot = y_dot_k - y_ref_dot_k;
    e_z_dot = z_dot_k - z_ref_dot_k;
    e_psi_dot = psi_dot_k - psi_ref_dot_k;
    
    % Initial u1 value for sliding mode coefficient calculation
    if k == 1
        u1_initial = m * g; % Initial estimate to balance gravity
    else
        u1_initial = U_history(1, k-1); % Use previous step's u1
    end
    
    % Update sliding mode coefficients according to Table 2
    % Check denominators to avoid numerical issues
    if abs(cos(phi_k) * cos(theta_k) * cos(psi_k)) > 1e-6
        c1 = 11 * m / (u1_initial * cos(phi_k) * cos(theta_k) * cos(psi_k));
        c2 = 6 * m / (u1_initial * cos(phi_k) * cos(theta_k) * cos(psi_k));
    else
        c1 = 11; % Default value when denominator is close to zero
        c2 = 6;
    end
    
    if abs(cos(psi_k)) > 1e-6
        c5 = -11 * m / (u1_initial * cos(psi_k));
        c6 = -6 * m / (u1_initial * cos(psi_k));
    else
        c5 = -11; % Default value when denominator is close to zero
        c6 = -6;
    end
    
    % Save sliding mode coefficient history
    c_history(k, :) = [c1, c2, c3, c4, c5, c6, c7, c8];
    
    % Define sliding surfaces
    s1 = c1 * e_x + c2 * e_x_dot;
    s2 = c5 * e_y + c6 * e_y_dot;
    s3 = e_z + cz * e_z_dot;
    s4 = c_psi * e_psi + c_psi * e_psi_dot;
    
    % Calculate virtual control inputs (for attitude reference generation)
    x_dot_dot_des = -(c1/c2) * e_x_dot - (epsilon_1/c2) * sat(s1, eta_1);
    y_dot_dot_des = -(c5/c6) * e_y_dot - (epsilon_2/c6) * sat(s2, eta_2);
    z_dot_dot_des = -cz * e_z_dot + epsilon_3 * sat(s3, eta_3) + g;
    psi_dot_dot_des = -c_psi * e_psi_dot - epsilon_4 * sat(s4, eta_4);
    
    % Save virtual control inputs for analysis
    virtual_control_history(:, k) = [x_dot_dot_des; y_dot_dot_des; z_dot_dot_des; psi_dot_dot_des];
    
    % Calculate desired attitude angles - Modified to ensure precision
    % Use small-angle approximation when accelerations are small to avoid numerical issues
    if abs(x_dot_dot_des) < 0.01 && abs(y_dot_dot_des) < 0.01
        phi_des = 0;
        theta_des = 0;
    else
        phi_des = asin((m * (x_dot_dot_des * sin(psi_k) - y_dot_dot_des * cos(psi_k))) / u1_initial);
        theta_des = asin((m * (x_dot_dot_des * cos(psi_k) + y_dot_dot_des * sin(psi_k))) / (u1_initial * cos(phi_des)));
    end
    
    % Limit desired attitude angles to safe range but allow more range
    phi_des = min(max(phi_des, -0.6), 0.6);  % About ±34 degrees
    theta_des = min(max(theta_des, -0.6), 0.6);  % About ±34 degrees
    
    % Save desired angles for analysis
    phi_des_history(k) = phi_des;
    theta_des_history(k) = theta_des;
    
    % Calculate attitude errors
    e_phi = phi_k - phi_des;
    e_theta = theta_k - theta_des;
    
    % Include reference angular velocities (derivative of desired angles)
    if k > 1
        phi_des_dot = (phi_des - phi_des_history(max(1, k-1))) / Ts;
        theta_des_dot = (theta_des - theta_des_history(max(1, k-1))) / Ts;
    else
        phi_des_dot = 0;
        theta_des_dot = 0;
    end
    
    e_phi_dot = phi_dot_k - phi_des_dot;
    e_theta_dot = theta_dot_k - theta_des_dot;
    
    % Define attitude control sliding surfaces
    s5 = c3 * e_phi + c4 * e_phi_dot;
    s6 = c7 * e_theta + c8 * e_theta_dot;
    
    % Calculate attitude control inputs
    phi_dot_dot_des = -(c3/c4) * e_phi_dot - (epsilon_3/c4) * sat(s5, eta_3);
    theta_dot_dot_des = -(c7/c8) * e_theta_dot - (epsilon_4/c8) * sat(s6, eta_4);
    
    % Calculate control inputs with improved model compensation
    u1 = m * (g + z_dot_dot_des + K3*z_dot_k/m) / (cos(phi_k) * cos(theta_k));
    
    % Torque control inputs - consider coupling terms and drag
    u2 = Ix * (phi_dot_dot_des + (theta_dot_k * psi_dot_k * (Iy - Iz) + Jr * theta_dot_k * Omega_r) / Ix + K4 * phi_dot_k / Ix);
    u3 = Iy * (theta_dot_dot_des + (phi_dot_k * psi_dot_k * (Iz - Ix) - Jr * phi_dot_k * Omega_r) / Iy + K5 * theta_dot_k / Iy);
    u4 = Iz * (psi_dot_dot_des + phi_dot_k * theta_dot_k * (Ix - Iy) / Iz + K6 * psi_dot_k / Iz);
    
    % Limit control inputs to safe range
    u1 = max(0.1, min(u1, 40));  % Thrust limit with minimum positive thrust
    u2 = max(-8, min(u2, 8));     % Roll torque limit (increased range)
    u3 = max(-8, min(u3, 8));     % Pitch torque limit (increased range)
    u4 = max(-8, min(u4, 8));     % Yaw torque limit (increased range)
    
    % Store control inputs
    U = [u1; u2; u3; u4];
    U_history(:, k) = U;
    
    % Simulate system response
    X = discreteC2D(X, U, Ts);
    
    % Store state
    X_history(:, k+1) = X;
end

% Transpose for plotting
X_history = X_history';
U_history = U_history';

% Create time vector for control inputs (one less point than tspan)
tspan_control = tspan(1:end-1);

% --- Plot Results ---
figure;

% Position tracking
subplot(3, 2, 1);
plot(tspan, X_history(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(tspan, x_ref, 'r--', 'LineWidth', 1.5);
legend('Actual', 'Reference');
title('Position x');
xlabel('Time (s)');
ylabel('x (m)');
grid on;

subplot(3, 2, 2);
plot(tspan, X_history(:, 3), 'b-', 'LineWidth', 2); hold on;
plot(tspan, y_ref, 'r--', 'LineWidth', 1.5);
legend('Actual', 'Reference');
title('Position y');
xlabel('Time (s)');
ylabel('y (m)');
grid on;

subplot(3, 2, 3);
plot(tspan, X_history(:, 5), 'b-', 'LineWidth', 2); hold on;
plot(tspan, z_ref, 'r--', 'LineWidth', 1.5);
legend('Actual', 'Reference');
title('Position z');
xlabel('Time (s)');
ylabel('z (m)');
grid on;

% Angles
subplot(3, 2, 4);
plot(tspan, X_history(:, 7), 'b-', 'LineWidth', 2); hold on;
plot(tspan_control, phi_des_history, 'r--', 'LineWidth', 1.5);
legend('Actual \phi', 'Desired \phi');
title('Roll Angle \phi');
xlabel('Time (s)');
ylabel('\phi (rad)');
grid on;

subplot(3, 2, 5);
plot(tspan, X_history(:, 9), 'b-', 'LineWidth', 2); hold on;
plot(tspan_control, theta_des_history, 'r--', 'LineWidth', 1.5);
legend('Actual \theta', 'Desired \theta');
title('Pitch Angle \theta');
xlabel('Time (s)');
ylabel('\theta (rad)');
grid on;

subplot(3, 2, 6);
plot(tspan, X_history(:, 11), 'b-', 'LineWidth', 2);
title('Yaw Angle \psi');
xlabel('Time (s)');
ylabel('\psi (rad)');
grid on;

% Control inputs
figure;
subplot(2, 2, 1);
plot(tspan_control, U_history(:, 1), 'g-', 'LineWidth', 2);
title('Thrust u1');
xlabel('Time (s)');
ylabel('Thrust (N)');
grid on;

subplot(2, 2, 2);
plot(tspan_control, U_history(:, 2), 'r-', 'LineWidth', 1.5);
title('Roll Torque u2');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;

subplot(2, 2, 3);
plot(tspan_control, U_history(:, 3), 'b-', 'LineWidth', 1.5);
title('Pitch Torque u3');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;

subplot(2, 2, 4);
plot(tspan_control, U_history(:, 4), 'g-', 'LineWidth', 1.5);
title('Yaw Torque u4');
xlabel('Time (s)');
ylabel('Torque (N·m)');
grid on;

% --- 3D Trajectory ---
figure;
plot3(X_history(:, 1), X_history(:, 3), X_history(:, 5), 'b-', 'LineWidth', 2);
hold on;
plot3(x_ref, y_ref, z_ref, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Trajectory');
legend('Actual', 'Reference');

% --- Virtual Control and Desired Angles Visualization ---
figure;
subplot(2, 2, 1);
plot(tspan_control, virtual_control_history(1,:), 'b-', 'LineWidth', 1.5);
title('x-axis Virtual Control (x_{dot\_dot\_des})');
xlabel('Time (s)');
ylabel('Acceleration (m/s²)');
grid on;

subplot(2, 2, 2);
plot(tspan_control, virtual_control_history(2,:), 'r-', 'LineWidth', 1.5);
title('y-axis Virtual Control (y_{dot\_dot\_des})');
xlabel('Time (s)');
ylabel('Acceleration (m/s²)');
grid on;

subplot(2, 2, 3);
plot(tspan_control, phi_des_history, 'b-', 'LineWidth', 1.5);
title('Desired Roll Angle \phi_{des}');
xlabel('Time (s)');
ylabel('\phi_{des} (rad)');
grid on;

subplot(2, 2, 4);
plot(tspan_control, theta_des_history, 'r-', 'LineWidth', 1.5);
title('Desired Pitch Angle \theta_{des}');
xlabel('Time (s)');
ylabel('\theta_{des} (rad)');
grid on;

% --- Sliding Mode Coefficient Visualization ---
figure;
subplot(2, 2, 1);
plot(tspan_control, c_history(:, 1), 'b-', tspan_control, c_history(:, 2), 'r--', 'LineWidth', 1.5);
title('Position x Sliding Mode Coefficients');
legend('c1', 'c2');
xlabel('Time (s)');
ylabel('Coefficient Value');
grid on;

subplot(2, 2, 2);
plot(tspan_control, c_history(:, 5), 'b-', tspan_control, c_history(:, 6), 'r--', 'LineWidth', 1.5);
title('Position y Sliding Mode Coefficients');
legend('c5', 'c6');
xlabel('Time (s)');
ylabel('Coefficient Value');
grid on;

subplot(2, 2, 3);
plot(tspan_control, c_history(:, 3), 'b-', tspan_control, c_history(:, 4), 'r--', 'LineWidth', 1.5);
title('Attitude \phi Sliding Mode Coefficients');
legend('c3', 'c4');
xlabel('Time (s)');
ylabel('Coefficient Value');
grid on;

subplot(2, 2, 4);
plot(tspan_control, c_history(:, 7), 'b-', tspan_control, c_history(:, 8), 'r--', 'LineWidth', 1.5);
title('Attitude \theta Sliding Mode Coefficients');
legend('c7', 'c8');
xlabel('Time (s)');
ylabel('Coefficient Value');
grid on;

% --- Helper Functions ---
% Saturation function (to reduce chattering)
function y = sat(s, phi)
    if abs(s) <= phi
        y = s/phi;
    else
        y = sign(s);
    end
end