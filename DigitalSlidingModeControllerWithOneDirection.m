clear; clc;

%% ---------------- 基本設定 ----------------
Ts   = 0.01;                  % 取樣時間
tend = 40;
t    = 0:Ts:tend;
N    = numel(t);

% 物理參數
m = 2; 
g = 9.8; 
K3 = 0.01;                    % z 方向線性阻尼
beta_hat_z = sqrt(0.2);       % 等效推力增益
beta_min_z = sqrt(0.2);       % 供你原式子使用（若不用可刪）
% 若只做 z 軸，先固定姿態 = 0（等效 cos(phi)*cos(theta)=1）
params.m = m; params.g = g; params.K3 = K3; params.beta_hat_z = beta_hat_z;
params.phi = 0; params.theta = 0;

%% ---------------- 參考軌跡 ----------------
z_ref     = ones(1, N);                     % 懸停 1 m
z_ref(round(N/2):end) = 2;                  % 中途切到 2 m
z_ref_dot = zeros(1, N);                    % 這裡取 0

%% ---------------- 陣列預配置 ----------------
z     = zeros(1, N);      % 位置
z_dot = zeros(1, N);      % 速度
U     = zeros(4, N);      % 控制輸入（只用 u1，其餘為 0）

%% ---------------- 初始條件 ----------------
z(1)     = 0;             
z(2)     = z(1);          % 二階差分需兩點
z_dot(1) = 0;
z_dot(2) = 0;

U(1,1) = m*g/beta_hat_z;             % 先用重力平衡
U(1,2) = m*g/beta_hat_z;

%% ---------------- SMC 參數（離散） ----------------
cz      = 1.0;            % 滑動面權重（s_z = cz*e_z + e_z_dot）
eta_z   = 5.0;            % 到達率（可調小一些減少抖動）
phi_bl  = 0.1;           % 邊界層厚度（0.02~0.1 可試）
u1_min  = 0.1; 
u1_max  = 100;

%% ---------------- 主迴圈（k = 2 ... N-1） ----------------
for k = 2:N-1
    % 1) ZOH：當步可用推力
    u1_now = U(1,k);

    % 2) 用當步狀態計算 a_z(k)
    a_z = accel_z( z(k), z_dot(k), u1_now, params );

    % 3) 二階差分更新 z
    z(k+1) = 2*z(k) - z(k-1) + Ts^2 * a_z;

    % 4) 速度離散更新（後向差分較穩）
    z_dot(k+1) = (z(k+1) - z(k))/Ts;

    % 5) SMC 誤差/滑動面（全用當步量）
    e_z     = z(k)     - z_ref(k);
    e_z_dot = z_dot(k) - z_ref_dot(k);
    s_z     = cz*e_z + e_z_dot;

    % 6) 邊界層 saturation 取代 sign
    sat_s = sat(s_z, phi_bl);

    % 7) 你的全離散控制律（把 sign 改成 sat）
    %    注意：分母包含姿態餘弦，這裡 cphi*cth = 1
    cphi = cos(params.phi); cth = cos(params.theta);
    denom = max(1e-6, beta_hat_z * cphi * cth);

    % 這條等同你之前的：u1_next = m*(g + (z(k+1)-z(k))*(K3/m - cz)/Ts - eta_z*sign(s_z)) / beta_hat_z
    u1_next = m*( g + (z(k+1)-z(k))*(K3/m - cz)/Ts - eta_z*sat_s ) / denom;

    % 8) 飽和並存到下一步（避免代數環）
    U(1, k+1) = min(max(u1_next, u1_min), u1_max);
end

%% ---------------- 繪圖 ----------------
figure; 
subplot(2,1,1);
plot(t, z, 'b', 'LineWidth', 1.6); hold on;
plot(t, z_ref, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('z (m)');
legend('z','z_{ref}'); title('Z position');

subplot(2,1,2);
plot(t, U(1,:), 'k', 'LineWidth', 1.0); grid on;
xlabel('Time (s)'); ylabel('u_1 (N)');
title('Thrust u_1');

%% ---------------- 輔助函數 ----------------
function a_z = accel_z(z, z_dot, u1, params)
    m = params.m; g = params.g; K3 = params.K3; beta_hat_z = params.beta_hat_z;
    cphi = cos(params.phi); cth = cos(params.theta);
    thrust = beta_hat_z * (u1/m) * cphi * cth;  % 推力等效加速度
    drag   = (K3/m) * z_dot;                    % 線性阻尼
    a_z = thrust - g - drag;
end

function y = sat(s, phi)
    % 邊界層 saturation：|s|<=phi 時線性，否則等於 sign(s)
    if abs(s) <= phi
        y = s/phi;
    else
        y = sign(s);
    end
end
