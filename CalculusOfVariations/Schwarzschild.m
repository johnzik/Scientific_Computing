clc;
clear;
tic

% --- Constants ---
Rs  = 1;

% --- Initial Conditions ---
r0 = 10;
phi0 = 0;
th0 = pi/2;

rFinal = 3;
phiFinal = pi/2;
thFinal = pi/2;

% dr0 = -4.5;
dphi0 = 30;
dth0 = 0;

options_fsolve = optimoptions('fsolve', 'Display', 'none');
dr0 = fsolve(@(test_dr0) get_r_error(test_dr0, r0, th0, phi0, phiFinal, rFinal, Rs, dth0), -4.5, options_fsolve);
fprintf('Ideal dr0 : %.5f\n', dr0);

% --- Hamiltonian Initial Conditions ---
qr0 = r0;
qphi0 = phi0;
qth0 = th0;

% Lagrangian
L = sqrt( dr0^2/((1 - Rs/r0)^2) + r0^2 * (dth0^2 + sin(th0)^2) / (1 - Rs/r0) );

pr0 = dr0 / ( (1 - Rs/r0)^2 * L );
pth0 = dth0 * r0^2 / ( (1 - Rs/r0) * L );

% --- Stricter tolerance ---
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

% --- Pack initial conditions ---
y0 = [qr0; qth0; pr0; pth0];

% --- Solve the 4 ODEs ---
phiSpan = linspace(phi0, phiFinal, 1000);
[phi_sol, Y_sol] = ode45(@(phi,y) F(phi, y, Rs), phiSpan, y0, options);

% --- Extract Solutions ---
rSol = Y_sol(:,1);
thSol = Y_sol(:,2);

% --- 3D Plot ---
x = rSol .* sin(thSol) .* cos(phi_sol);
y = rSol .* sin(thSol) .* sin(phi_sol);
z = rSol .* cos(thSol);
% [xcart, ycart] = pol2cart(thSol, rSol);

figure;
hold on;
grid on;
axis equal;

% Plot Black Hole
% th_circle = linspace(0, 2*pi, 100);
% x_bh = Rs * cos(th_circle);
% y_bh = Rs * sin(th_circle);
% fill(x_bh, y_bh, 'k');
% plot3(0, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
[Xbh, Ybh, Zbh] = sphere(50); 
surf(Xbh, Ybh, Zbh, 'FaceColor', 'k', 'EdgeColor', 'none');


% Takeoff and Landing positions
plot(x(1), y(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5); % Earth
plot(x(end), y(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); % Miller Black Hole

title('Trajectory of Endurance Spacecraft');

% Plot Null Geodesic
plot(x, y, 'b', 'LineWidth', 2);
legend('Event Horizon (Rs=1)', 'Takeoff (10, 0)', 'Landing (0, 3)', 'Null Geodesic');
view(3);
% comet(x,y);

disp('Done');
toc

function dy = F(t, y, Rs)
    % Unpack y0
    r = y(1);
    th = y(2);
    pr = y(3);
    pth = y(4);
    
    f = 1 - Rs/r;
    
    dy = zeros(4,1);                                                                    % Initialize output vector
    dy(1) = r * sin(th) * pr * f / sqrt( 1/f - pr^2*f - pth^2/r^2 );                    % qr'
    dy(2) = pth * sin(th) / ( r * sqrt( 1/f - pr^2*f - pth^2/r^2 ) );                   % qth'
    dy(3) = sin(th)*(pr^2*(Rs/r - 1) - pth^2/r^2 - 1/(Rs/r - 1))^(1/2) - ...
            (r*sin(th)*((Rs*pr^2)/r^2 - (2*pth^2)/r^3 + Rs/(r^2*(Rs/r - 1)^2)))/ ...
            (2*(pr^2*(Rs/r - 1) - pth^2/r^2 - 1/(Rs/r - 1))^(1/2));                     % pr'
    dy(4) = r*cos(th)*(pr^2*(Rs/r - 1) - pth^2/r^2 - 1/(Rs/r - 1))^(1/2);               % pth'
end

function err = get_r_error(dr0_test, r0, th0, phi0, phiFinal, rFinal, Rs, dth0)
    L = sqrt( dr0_test^2/((1 - Rs/r0)^2) + r0^2 * (dth0^2 + sin(th0)^2) / (1 - Rs/r0) );
    pr0 = dr0_test / ( (1 - Rs/r0)^2 * L );
    pth0 = dth0 * r0^2 / ( (1 - Rs/r0) * L );
    y0 = [r0; th0; pr0; pth0];
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    
    [~, Y_test] = ode45(@(t,y) F(t, y, Rs), [phi0, phiFinal], y0, options);
    
    err = Y_test(end,1) - rFinal;
end

% syms r Rs pr pth th
% 
% H = -r*sin(th)*sqrt( 1/(1-(Rs/r)) - pr^2*(1-(Rs/r)) - pth^2/r^2 );
% 
% dpr = -diff(H,r);
% dpth = -diff(H,th);