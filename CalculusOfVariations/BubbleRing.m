%% A) Find y(x)
clc;
clear;

% Initialize the parameters
R1 = 1;
R2 = 0.8;

L = 1;                  % Distance between the two rings     
S = pi*R1^2 + pi*R2^2;  % Surface of the 2 rings combined
step = 0.1;             % Step that increases the distance between the rings

xmesh = linspace(0, L, 100);
solinit = bvpinit(xmesh, @guess);

% Solve the boundary value problem
sol = bvp4c(@odefun, @bcfun, solinit);

% Plot the solution
x = sol.x;
y = sol.y(1, :);
yp = sol.y(2, :);

figure;
plot(x, y, 'b-'); grid;
hold on;
plot(x, -y, 'b-');

plot([0 0], [R1 -R1], 'k-'); 
plot([1 1], [R2 -R2], 'k-');
title('2D Line Solution of the min y(x)')

% --- 3D plot ---
theta = linspace(0, 2*pi, 100);

[THETA, x_3D] = meshgrid(theta, x);

y_radius = repmat(y', 1, length(theta));

y_3D = y_radius .* cos(THETA);
z_3D = y_radius .* sin(THETA);

figure;
surf(x_3D, y_3D, z_3D, 'EdgeColor', 'none', 'FaceAlpha', 0.6); title('3D Solution Bubble');
hold on;

% Plot the rings
plot3(zeros(1,100), R1*cos(theta), R1*sin(theta), 'k-', 'LineWidth', 3);
plot3(L*ones(1,100), R2*cos(theta), R2*sin(theta), 'k-', 'LineWidth', 3); grid on;

% --- Analytic Solution ---

syms y(x)
DE = diff(y,x,2) == ( 1 + diff(y,x)^2 ) / y;

sol = dsolve(DE);

disp('The analytic equation is y(x):');
disp(sol);

%% B) Which is the Lcrit where the bubble bursts?
clear;

% Initialize the radiuses
R1 = 1;
R2 = 0.8;

L = 0.1;                  % Distance between the two rings     
S = pi*R1^2 + pi*R2^2;  % Surface of the 2 rings combined
step = 0.01;             % Step that increases the distance between the rings

xmesh = linspace(0, L, 100);
solinit = bvpinit(xmesh, @guess);

% Solve the boundary value problem
sol = bvp4c(@odefun, @bcfun, solinit);

% Plot the solution
x = sol.x;
y = sol.y(1, :);
yp = sol.y(2, :);

integrand = y .* sqrt( 1 + yp.^2);
surf_3D = 2*pi * trapz(x, integrand);

while true
    % Break if the area if surf_3D is larger than S
    if surf_3D >= S
        disp(['The bubble will burst since surf_3D = ', num2str(surf_3D), ' > S = ', num2str(S)]);
        disp(['Lcrit = ', num2str(L)]);
        break;
    end
    
    L = L + step;
    
    xmesh = linspace(0, L, 100);
    solinit = bvpinit(xmesh, @guess);

    % Solve the boundary value problem
    sol = bvp4c(@odefun, @bcfun, solinit);

    x = sol.x;
    y = sol.y(1, :);
    yp = sol.y(2, :);    
    
    % Compute the new area
    integrand = y .* sqrt( 1 + yp.^2);
    surf_3D = 2*pi * trapz(x, integrand);   
end

function dydx = odefun(x, y)
    dydx = zeros(2,1);
    dydx(1) = y(2);                     % y1' = y2
    dydx(2) = (1 + y(2)^2) / y(1);      % y2' = (1 + y2^2) / y1
end

function res = bcfun(ya, yb)
    res = [ ya(1) - 1;                  % y(0) = 1
            yb(1) - 0.8];               % y(1) = 0.8 
end

function yinit = guess(x)
    yinit = [1; 0]; 
end