clc;
clear; 

% Constants
A0 = 0.6;  % A0 = x m^2 is the cross-section area under y(x)

% Initial guess
x0 = [-0.1; -0.1];

x = fsolve(@(x) F(x, A0), x0);

% Display results
disp(['Lambda: ' num2str(x(1))])
disp(['c_2: ' num2str(x(2))])

% Check results
disp(F(real(x), A0));

% Create plot
x = real(x);
t = 0:0.001:1;
y = (-1/x(1)) * sqrt( 1 - (x(1)*t).^2 ) + x(2);
plot(t,y, 'LineWidth', 3);
xlabel('x');
ylabel('y');
title(['Optimal Arc for Venturi Vent with Cross-Sectional Area A0: ', num2str(A0), 'm^2']);
hold on;
grid on;
axis equal;
ylim([0, inf]);

% plot([0, 1], [0.06, 0.06]);
hold off;

% x1 = lamda, x2 = c2
function y = F(x, A0)
    y(1) = ( -1/x(1) ) * sqrt( 1 - x(1)^2 ) + x(2) - 0.06;
    y(2) = x(2) - (asin(x(1))/2 + (x(1)*(1 - x(1)^2)^(1/2))/2)/x(1)^2 - A0;
end