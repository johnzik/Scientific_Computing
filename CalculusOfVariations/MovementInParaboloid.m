clc;
clear;

% --- Constants ---
m = 10;
g = 9.8;
Tmax = 10;

% --- Initial conditions ---
r0 = 2;
phi0 = pi/2;
z0 = 4;
dr0 = 2;
dphi0 = 4.43; % dphi0 = 4.43 & dr0 = 0 for circle around z-axis
dz0 = 0;

% --- Hamiltonian initial conditions ---
qr_0 = r0;
qphi_0 = phi0;
qz_0 = z0;

pr_0 = m * dr0 * (1+4*r0^2);
pphi_0 = m * r0^2 * dphi0;
pz_0 = 0;

% --- Stricter tolerance ---
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

% --- Pack initial conditions ---
y0 = [qr_0; qphi_0; pr_0; pphi_0];

% --- Solve the 6 ODEs ---
tspan = linspace(0, Tmax, 1000);
[t, Y] = ode45(@(t,y) F(t, y, m, g), tspan, y0, options);

% --- Extract the data ---
r_out   = Y(:, 1);
phi_out = Y(:, 2);

% --- 3D Plot ---
x_sol = r_out .* cos(phi_out);
y_sol = r_out .* sin(phi_out);
z_sol = r_out.^2;

hold on;
% plot3(x_sol, y_sol, z_sol, 'r', 'LineWidth', 2);
grid on;

% Plot a Paraboloid  
x = linspace(-3,3,100);
y = linspace(-3,3,100);
[parabX, parabY] = meshgrid(x,y);
Fparab = parabX.^2 + parabY.^2;
surf(parabX,parabY,Fparab, 'FaceAlpha', '0.5', 'EdgeColor', 'none');
title('Trajectory of a Mass Inside a Paraboloid');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(3);

info_text = sprintf('Initial Velocities:\nv_r = %.2f m/s\nv_ö = %.2f rad/s\nv_z = %.2f m/s', dr0, dphi0, dz0);

annotation('textbox', [0.15 0.75 0.2 0.15], 'String', info_text, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'FaceAlpha', 0.8, ...
    'EdgeColor', 'black', ...
    'FontWeight', 'bold');

plot3(0, 2, 4, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

% Plot the trajectory
comet3(x_sol, y_sol, z_sol);
grid on;

function dy = F(t, y0, m, g)
    % Unpack y0
    r    = y0(1);
    phi   = y0(2);
    pr   = y0(3);
    pphi  = y0(4);
    
    dy = zeros(4,1);                                                % Initialize output vector
    dy(1) = pr / (m * (1+4*r^2));                                   % q_r' 
    dy(2) = pphi / (m*r^2);                                         % q_phi'
    dy(3) = (4*r*pr^2) / (m * (1+4*r^2)^2) + pphi^2 / (m*r^3) - ...
             2*m*g*r;                                               % p_r'
    dy(4) = 0;                                                      % p_phi'
end