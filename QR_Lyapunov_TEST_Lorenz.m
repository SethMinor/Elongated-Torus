% Testing Lyapunov exponent code on Lorenz system
clear, clc, clf;

% Font size, for plotting
fs = 10;
%'Interpreter','latex','FontSize', fs

% Lorenz system, Lyapunov spectrum
% lambda_1 ~ 0.9056
% lambda_2 = 0
% lambda_3 ~ -14.572

% Parameters
p = 10;
r = 28;
b = 8/3;

% Integrate the equations of motion
% Set total time and tolerances
t0 = 0;
tf = 500;
timespan = [t0, tf];
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Numerical integration using ode45
y0 = [0.1, 0.1, 0.1];
[t,y] = ode45(@(t,y) LYAP_RHS(t,y,p,r,b),timespan, y0, options);

figure (1)
plot3(y(:,1), y(:,2), y(:,3))
grid on

%% Compute Lyapunov exponents
% Roll through 'y_n' solution and store {r^n} values
skip_every = 1;

% Define the RHS
F =@(Y) LYAP_RHS(0,[Y(1), Y(2), Y(3)], p, r, b);

% Compute variables for dQ/dt = QS
Q_n = eye(3); % Initialize Q_0
J_0 = myjacobian(F, y(1,:)); % Initialize Jacobian matrix

% Define S matrix
% triu + tril commands
temp = (Q_n)'*J_0*Q_n;
S_0 = triu(-temp') + tril(temp);

% Solve dQ/dt = QS
Q_list = zeros(3,3,length(y)-1);
Q_list(:,:,1) = Q_n;

% RK-4 by hand
for n = 2:length(y)-1
    dt = t(n) - t(n-1);

    k1 = RHS_Q(Q_n,F,y,n-1,0);
    k2 = RHS_Q(Q_n + dt*k1/2,F,y,n-1,1);
    k3 = RHS_Q(Q_n + dt*k2/2,F,y,n-1,1);
    k4 = RHS_Q(Q_n + dt*k3,F,y,n,0);
    
    % Solve
    Q_new = Q_n + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    Q_n = Q_new;
    
    % Store Q_n solution
    Q_list(:,:,n) = Q_n;
end

% Solve {rho_ii} ODE
rho_n = ones(1,3); % Initial condition
rho_list = zeros(length(y)-1,3);
rho_list(1,:) = rho_n;

% RK-4 by hand (here we go again)
for n = 2:length(y)-1
    dt = t(n) - t(n-1);
    
    % Find Q_n and Q_{n+1}
    Q_n = Q_list(:,:,n-1);
    Q_nplus = Q_list(:,:,n);

    % Ghetto RK-4 (?)
    k1 = RHS_rho(Q_n,F,y,n-1,0);
    k2 = RHS_rho((Q_n + Q_nplus)/2,F,y,n-1,1);
    k3 = k2;
    k4 = RHS_rho(Q_nplus,F,y,n,0);
    
    % Solve
    rho_new = rho_n + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    rho_n = rho_new;
    
    % Store Q_n solution
    rho_list(n,:) = rho_n;
end

% Compute the full spectrum
N_lyapunov = 200; % number to back-average
Lyapunov = zeros(1,3);
tnew = t(2:end); % To make things easier
for i = 1:3
    Lyapunov(i) = mean(rho_list(end-N_lyapunov:end,i)./tnew(end-N_lyapunov:end));
end

% Accepted values
L1 = 0.9056;
L2 = 0;
L3 = -14.572;

%% Plot the spectrum
figure (2)
sgtitle('Lyapunov Spectrum (Lorenz)','Interpreter','latex','FontSize',fs+2)

subplot(3,1,1)
plot(tnew, rho_list(:,1)./tnew)
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
plot(tnew, rho_list(:,2)./tnew)
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
plot(tnew, rho_list(:,3)./tnew)
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
function dydt = LYAP_RHS(~, y, p, r, b)
  dydt(1) = p*(y(2) - y(1));
  dydt(2) = (r - y(3))*y(1) - y(2);
  dydt(3) = -b*y(3) + y(1)*y(2);

  dydt = dydt';
end

% Numerical Jacobian matrix
% See (https://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf)
function J = myjacobian(func,x)
    % Set numerical derivative parameters
    N = length(x);
    %F_at_x = feval(func,x);
    epsilon = 1E-8;

    % Compute numerical derivative
    xperturb = x;
    xperturb2 = x; 
    xperturb_minus = x;
    xperturb_minus2 = x;
    %xperturb = x + epsilon;
    J = zeros(N);
    for i = 1:N
        xperturb(i) = xperturb(i) + epsilon;
        xperturb2(i) = xperturb2(i) + 2*epsilon;

        xperturb_minus(i) = xperturb_minus(i) - epsilon;
        xperturb_minus2(i) = xperturb_minus2(i) - 2*epsilon;

        %J(:,i) = (feval(func,xperturb) - F_at_x)/epsilon;
        %J(:,i) = (feval(func,xperturb) - feval(func,xperturb_minus))/(2*epsilon);
        J(:,i) = (-feval(func,xperturb2) + 8*feval(func,xperturb) ...
            - 8*feval(func,xperturb_minus) + feval(func,xperturb_minus2))/(12*epsilon);

        xperturb(i) = x(i);
        xperturb2(i) = x(i);
        xperturb_minus(i) = x(i);
        xperturb_minus2(i) = x(i);
    end
end

% Modified Gram-Schmidt QR, for Lyapunov exponents
% https://www.mathworks.com/matlabcentral/fileexchange/55881-gram-schmidt-orthogonalization
function [Q, R] = Gram_Schmidt_QR(X)
    % Modified Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm) 
    % which produces the same result as [Q,R]=qr(X,0)
    % Written by Mo Chen (sth4nth@gmail.com).
    [d,n] = size(X);
    m = min(d,n);
    R = zeros(m,n);
    Q = zeros(d,m);
    for i = 1:m
        v = X(:,i);
        for j = 1:i-1
            R(j,i) = Q(:,j)'*v;
            v = v-R(j,i)*Q(:,j);
        end
        R(i,i) = norm(v);
        Q(:,i) = v/R(i,i);
    end
    R(:,m+1:n) = Q'*X(:,m+1:n);
end

% dQ/dt = QS ODE for Lyapunov Spectrum
function dQdt = RHS_Q(Q,F,y,n,midpoint_bool)
  % Jacobian at y_n
  J_n = myjacobian(F, y(n,:)); % Initialize Jacobian matrix

  % Jacobian at y_{n+1}
  J_nplus = myjacobian(F, y(n+1,:));

  % Decide on a Jacobian averaging
  if midpoint_bool == 0
    temp = (Q')*J_n*Q;
  elseif midpoint_bool == 1
    Jtemp = (J_n + J_nplus)/2;
    temp = (Q')*Jtemp*Q;
  end

  % Define the S matrix
  S = triu(-temp') + tril(temp);

  % Return RHS 
  dQdt = Q*S;
end

% ODE for Lyapunov Spectrum, d(rho_ii)/dt = (Q'JQ)_ii
function dpdt = RHS_rho(Q,F,y,n,midpoint_bool)
  % Jacobian at y_n
  J_n = myjacobian(F, y(n,:)); % Initialize Jacobian matrix

  % Jacobian at y_{n+1}
  J_nplus = myjacobian(F, y(n+1,:));

  % Decide on a Jacobian averaging
  if midpoint_bool == 0
    temp = (Q')*J_n*Q;
  elseif midpoint_bool == 1
    Jtemp = (J_n + J_nplus)/2;
    temp = (Q')*Jtemp*Q;
  end

  % Return RHS 
  dpdt = diag(temp)';
end
