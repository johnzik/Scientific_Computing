clc;
clear;

% --- Constants ---
m = 1;
l0 = 1;
g = 9.8;
k = 50;
Tmax = 10;

% --- Initial conditions ---
r0 = 2;
th0 = pi/6;
phi0= pi/4;
dr0 = 1;
dth0 = 2;
dphi0 = 2;

% --- Hamiltonian initial conditions ---
qr_0 = r0;
qth_0 = th0;
qphi_0 = phi0;
pr_0 = m*dr0;
pth_0 = m*(r0^2)*dth0;
pphi_0 = m * (r0^2) * (sin(th0)^2) * dphi0;

% --- Stricter tolerance ---
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

% --- Pack initial conditions ---
y0 = [qr_0; qth_0; qphi_0; pr_0; pth_0; pphi_0];

% --- Solve the 6 ODEs ---
[t, Y] = ode45(@(t,y) F(t, y, m, g, l0, k), [0 Tmax], y0, options);

% --- Extract the data ---
r_out   = Y(:, 1);
th_out  = Y(:, 2);
phi_out = Y(:, 3);

% --- 3D Plot ---
x = r_out .* sin(th_out) .* cos(phi_out);
y = r_out .* sin(th_out) .* sin(phi_out);
z = -r_out .* cos(th_out); 

figure;
plot3(x, y, z, 'b', 'LineWidth', 1.5);
grid on;
title('Trajectory of Sping Pendulum');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

hold on;
plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend('Mass Trajectory', 'Pivot Point (V=0)');

figure;
comet3(x, y, z);
grid on;
title('Trajectory of Sping Pendulum');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

function dy = F(t, y0, m, g, l0, k)
    % Unpack y0
    r    = y0(1);
    th   = y0(2);
    phi  = y0(3);
    pr   = y0(4);
    pth  = y0(5);
    pphi = y0(6);
    
    % Initialize output vector
    dy = zeros(6,1);
    dy(1) = pr / m;                                             % q_r' 
    dy(2) = pth / (m*r^2);                                      % q_th'
    dy(3) = pphi / (m*r^2*sin(th)^2);                           % q_phi'
    dy(4) = pth^2 / (m*r^3) + pphi^2 / (m*r^3*sin(th)^2) + ...
            m*g*cos(th) - k*(r-l0);                             % p_r'
    dy(5) = (pphi^2 / (m*r^2)) * csc(th)^2 * cot(th) - ...
            m*g*r*sin(th);                                      % p_th'
    dy(6) = 0;                                                  % p_phi'
end