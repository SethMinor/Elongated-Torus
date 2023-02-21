% Testing Lyapunov exponent code on Lorenz system
clear, clc, clf;

% Font size, for plotting
fs = 10;
%'Interpreter','latex','FontSize', fs

% Lorenz system, Lyapunov spectrum
% Accepted values
L1 = 0.9056;
L2 = 0;
L3 = -14.572;

% Parameters
p = 10;
r = 28;
b = 8/3;

% Integrate the equations of motion
% Set total time and tolerances
t0 = 0;
tf = 100;
timespan = [t0, tf];
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Define the RHS of the Lorenz system
% For computing the numerical Jacobian
F =@(y) LORENZ_RHS(0,[y(1), y(2), y(3)], p, r, b);

% Numerical integration using ode45
% RHS = [dydt (N); dQdt (N^2); dpdt (N)]
y0 = [0.1, 0.1, 0.1]'; % Lorenz ICs
Q0 = eye(3); % Lyapunov, Q ICs
p0 = ones(3,1); % Lyapunov, rho_i ICs
Y0 = [y0; reshape(Q0,[length(Q0)^2,1]); p0];

[t,y] = ode45(@(t,y) LYAPUNOV_RHS(t,y, p,r,b, F), timespan, Y0, options);

% Plot the attractor
figure (1)
plot3(y(:,1), y(:,2), y(:,3))
grid on

%% Compute Lyapunov exponents
% Roll through 'y_n' solution and store {rho_i} values
% [rho_1, ..., rho_4] list = y(:, 1 + N + N^2 : end)
rho_list = y(:, 1+length(y0)+length(y0)^2 : end);

% Compute the full spectrum
N_lyapunov = 200; % number to back-average
Lyapunov = zeros(1,3);
for i = 1:3
    Lyapunov(i) = mean(rho_list(end-N_lyapunov:end,i)./t(end-N_lyapunov:end));
end

%% Plot the spectrum
plot_every = 10; % Speed-up for plotting

figure (2)
sgtitle('Lyapunov Spectrum (Lorenz)','Interpreter','latex','FontSize',fs+2)

subplot(3,1,1)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,1)./t(1:plot_every:end))
hold on
yline(Lyapunov(1),'--k')
yline(L1,'-r','LineWidth',0.3)
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{1}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_1=$ "+Lyapunov(1),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Accepted value','Interpreter','latex')
ylim([0.5, 1.5])

subplot(3,1,2)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,2)./t(1:plot_every:end))
hold on
yline(Lyapunov(2),'--k')
yline(L2,'-r','LineWidth',0.3)
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{2}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_2=$ "+Lyapunov(2),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Accepted value','Interpreter','latex')
ylim([-0.5, 0.5])

subplot(3,1,3)
plot(t(1:plot_every:end), rho_list(1:plot_every:end,3)./t(1:plot_every:end))
hold on
yline(Lyapunov(3),'--k')
yline(L3,'-r','LineWidth',0.3)
hold off
grid on
xlabel('$t$','Interpreter','latex','FontSize',fs)
ylabel('$\rho_{3}/t$','Interpreter','latex','FontSize',fs)
title("$\lambda_3=$ "+Lyapunov(3),'Interpreter','latex','FontSize',fs)
legend('Numerical','Averaged','Accepted value','Interpreter','latex')
ylim([-15, -14])

%% Function definitions
% RHS of Lorenz system
function dydt = LORENZ_RHS(~, y, p, r, b)
  % Lorenz system eq's of motion
  dydt(1) = p*(y(2) - y(1));
  dydt(2) = (r - y(3))*y(1) - y(2);
  dydt(3) = -b*y(3) + y(1)*y(2);
  dydt = dydt';
end

% RHS for entire Lyapunov spectrum
function dYdt = LYAPUNOV_RHS(~,y, p,r,b, F)
  % Lorenz system eq's of motion
  dydt(1) = p*(y(2) - y(1));
  dydt(2) = (r - y(3))*y(1) - y(2);
  dydt(3) = -b*y(3) + y(1)*y(2);
  dydt = dydt';

  % Compute Jacobian for later use
  J = myjacobian(F, y(1:length(dydt))); % Jacobian at y

  % Integrate dQ/dt = QS ODE
  % Reshape Q into matrix, Q = y(1+N : N+N^2)
  Q = y(length(dydt)+1 : length(dydt)+length(dydt)^2);
  Q = reshape(Q, [length(dydt), length(dydt)]);
  % Compute S
  temp = (Q')*J*Q;
  S = triu(-temp') + tril(temp);
  % Reshape RHS as vector
  dQdt = Q*S;
  dQdt = reshape(dQdt, [length(dydt)^2, 1]);

  % Integrate d(rho_i)/dt = (Q'JQ)_ii ODE
  % rho_i = y(1 + N + N^2 : end)
  % Return RHS 
  dpdt = diag(temp);

  % Return the combined RHS
  dYdt = [dydt; dQdt; dpdt];
end

% Numerical Jacobian matrix
% See (https://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf)
function J = myjacobian(func,x)
    % Set numerical derivative parameters
    N = length(x);

    %F_at_x = feval(func,x);
    epsilon = 1E-8; % Step size

    % Compute numerical derivative
    xperturb = x;
    xperturb2 = x; 
    xperturb_minus = x;
    xperturb_minus2 = x;
    %xperturb = x + epsilon;
    J = zeros(N);

    for i = 1:N
        % Perturbation from current state
        xperturb(i) = xperturb(i) + epsilon; % x + h
        xperturb2(i) = xperturb2(i) + 2*epsilon; % x + 2h

        xperturb_minus(i) = xperturb_minus(i) - epsilon; % x - h
        xperturb_minus2(i) = xperturb_minus2(i) - 2*epsilon; % x - 2h
        
        % First order scheme
        %J(:,i) = (feval(func,xperturb) - F_at_x)/epsilon;

        % Second order scheme
        %J(:,i) = (feval(func,xperturb) - feval(func,xperturb_minus))/(2*epsilon);

        % Third order scheme
        J(:,i) = (-feval(func,xperturb2) + 8*feval(func,xperturb) ...
            - 8*feval(func,xperturb_minus) + feval(func,xperturb_minus2))/(12*epsilon);
        
        % Reset loop variables
        xperturb(i) = x(i);
        xperturb2(i) = x(i);
        xperturb_minus(i) = x(i);
        xperturb_minus2(i) = x(i);
    end
end
