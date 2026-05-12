clc;
clear;
tic
% --- Initial Conditions & Z Limit ---
% y0 = [r(0), r'(0)];
y0 = [0.05; 0.5;];
L = 1;
zspan = [0, L];
% a = 50;
lamda0 = 600;   % 600nm
lamda_r = 633;  % 633nm for HeNe
lamda_g = 514;  % 514nm for Argon
a = [20*( 1 + (lamda0/lamda_r)^2 ), 20*( 1 + (lamda0/lamda_g)^2 )];
color = ['r', 'g'];

figure;
hold on;
for i = 1:2
    
    % --- Set up the 2 ODEs ---
    % Set y1 = r, y2 = r'
    odesys = @(z, y) [
        y(2);                       % y1' = r' = y2
        -2*a(i)*y(1)*(1+y(2)^2);    % y2' = r''
    ];
    
    % --- Solve with ode45 ---
    [z_sol, y_sol] = ode45(odesys, zspan, y0);
    
    % --- 3D Plot ---
    X_val = y_sol(:,1);
    
    Y_val = zeros(size(z_sol)); 
    
    plot3(z_sol, X_val, Y_val, color(i), 'LineWidth', 1.5);
end

[Xc, Yc, Zc] = cylinder(0.1, 50);
Zc = Zc * L;
surf(Zc, Xc, Yc, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off', 'FaceColor', 'c');

grid on;
xlabel('Z (m) - Cylindrical Waveguide Length');
ylabel('X (m) - Radius r');
zlabel('Y (m)');
title(['Meridional Rays in GRIN Waveguide with a_r: ', num2str(a(1)), ' and a_g', num2str(a(2))]);
legend('HeNe (\lambda_r = 633 nm)', 'Argon (\lambda_g = 514 nm)', 'Location', 'best');
axis equal;
zlim([-0.1 0.1]);
ylim([-0.1 0.1]);
view(3);

disp('Done!');
toc