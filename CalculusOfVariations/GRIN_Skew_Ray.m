clc;
clear;
tic

% --- Initial Conditions & Z Limit ---
% y0 = [r(0), r'(0), th(0), th'(0)];
y0 = [0.05; 0.5; 0; 6];
L = 1;
zspan = [0, L];
% a = 50;
lamda0 = 600;   % 600nm
lamda_r = 633;  % 633nm for HeNe
lamda_g = 514;  % 514nm for Argon
a = [20*( 1 + (lamda0/lamda_r)^2 ), 20*( 1 + (lamda0/lamda_g)^2 )];
color = ['r', 'g'];

syms r(z) th(z)
dr = diff(r,z);
dth = diff(th,z);

figure;
hold on;
for i = 1:2
    f = exp(-a(i)*r^2) * sqrt(1 + dr^2 + r^2 * dth^2);

    % Auto EL calc: df/dr - d/dz(df/dr')
    EL_r = functionalDerivative(f, r) == 0;

    % Auto EL calc: df/dth - d/dz(df/dth')
    EL_th = functionalDerivative(f, th) == 0;

    % --- Solve for the 2nd order variables ---
    syms ddr_dummy ddth_dummy

    EL_r_sub = subs(EL_r, diff(r, z, 2), ddr_dummy);
    EL_th_sub = subs(EL_th, diff(th, z, 2), ddth_dummy);

    [sol_ddr, sol_ddth] = solve([EL_r_sub, EL_th_sub], [ddr_dummy, ddth_dummy]);

    % --- Swap r(z), diff() with sym variables ---
    syms r_val th_val dr_val dth_val

    eq_r = subs(sol_ddr, [r, dr, th, dth], [r_val, dr_val, th_val, dth_val]);
    eq_th = subs(sol_ddth, [r, dr, th, dth], [r_val, dr_val, th_val, dth_val]);

    % --- Convert to anonymous functions ---
    fddr = matlabFunction(eq_r, 'Vars', [r_val, dr_val, th_val, dth_val]);
    fddth = matlabFunction(eq_th, 'Vars', [r_val, dr_val, th_val, dth_val]);

    % --- Set up the 4 ODEs ---
    % Set y1 = r, y2 = r', y3 = th, y4 = th'
    ode_system = @(z, y) [
        y(2);                           % y1' = r' = y2
        fddr(y(1), y(2), y(3), y(4));   % y2' = r'' = fddr
        y(4);                           % y3' = th' = y4
        fddth(y(1), y(2), y(3), y(4));  % y4' = th'' = fddth
    ];

    % --- Solve with ode45 ---
    [z_sol, y_sol] = ode45(ode_system, zspan, y0);

    % --- 3D Plot ---
    r_sol = y_sol(:, 1);
    th_sol = y_sol(:, 3);

    % Convert Cylindrical Coords to Cartesian Coords
    X = r_sol .* cos(th_sol);
    Y_cartesian = r_sol .* sin(th_sol);

    plot3(z_sol, X, Y_cartesian, 'LineWidth', 2, 'Color', color(i));
end

[Xc, Yc, Zc] = cylinder(0.1, 50);
Zc = Zc * 1;
surf(Zc, Xc, Yc, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'c');

xlabel('Z (m) - Cylindrical Waveguide Length');
ylabel('X (m) - Radius r');
zlabel('Y (m)');
title(['Skew Rays in GRIN Waveguide with a_r: ', num2str(a(1)), ' and a_g', num2str(a(2))]);
legend('HeNe (\lambda_r = 633 nm)', 'Argon (\lambda_g = 514 nm)', 'Location', 'best');
grid on;
axis equal;
zlim([-0.1 0.1]);
ylim([-0.1 0.1]);
view(3);

disp('Done!');
toc