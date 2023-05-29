%% Branching for MLE plot, v2 (speed-up)
% This script computes a branch of points under torus deformation.
% Also allows for computation of MLE using branch points as ICs.
clear, clc, clf;

% Font size, for plotting
fs = 18;

% Compute Lyapunov exponents of solution branch?
Lyap_bool = true;

% Initial parameters
a = 13; % Current semi-major axis
a0 = a; % Initial semi-major axis
R = 11;
r = 3;
c = sqrt(R^2 - r^2);
myalpha = r/R;

% Energy level to check
E0 = -9.87844;

% Initial seeds for [v1,v2]
v_seed = [-10, -9.8];

% Set u2 = constant = u1
u2constant = -30;

%% Select a branch of ICs for testing
% Continuation for-loop
da = 0.25; % delta-a step (semi-major axis step)
counter = 1;

% Variables used in the loop
delta = 1E-6; % MLE computation
ExpFac = 1E2; % MLE computation
cap = 20; % Jacobi-theta truncation
phi_0 = 0; % Isothermal boundary condition
theta_span = [0, pi];
q = [1, -1];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

for delta_a = 0:da:1 % ('a' will go from a0 to a0-2)
    % New semi-major axis a
    a = a0 - delta_a;

    % Initial seed
    if counter == 1
        V = v_seed;
    else
        V = v_star;
    end
    
    % --- Solving for isothermal coordinates ---
    % Numerical integration using ode45
    options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);
    [theta_raw,phi_raw] = ode15s(@(theta,phi) odefcn(theta,phi,a,R,r), theta_span, phi_0, options);
    
    % Give these hosses odd symmetry 
    theta_raw = [-flip(theta_raw(2:end)); theta_raw];
    phi_raw = [-flip(phi_raw(2:end)); phi_raw];
    
    % Define gl and gr
    gl = min(imag(phi_raw));
    gr = max(imag(phi_raw));
    
    % Defining the (u,v) to (phi,theta) map
    ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
    phi =@(u) u/c;
    theta =@(v) ginv(-UVwrap(v, [-c*gr, c*gr])/c);

    % Define the derivative of g_inverse
    D2 = differentiate(ginv, imag(phi_raw));
    Dginv = fit(imag(phi_raw), D2, 'cubicinterp');
    
    % --- Other parameters ---
    % Gamma quantity
    gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
    
    % Jacobi theta function parameters
    p = exp(-gr); % nome (for periodicity)

    % Define implicit function f = H(w1,w2)-E0
    f =@(u1,v1,v2,E0) log(abs( jacobitheta1(((u1 - u2constant) + 1i*(v1-v2))./(2*c),p,cap) )) + ...
        log(abs( jacobitheta1(-((u1 - u2constant) + 1i*(v1-v2))./(2*c),p,cap) )) + ...
        log(gamma(u1./c, theta(v1)).*gamma(u2constant./c, theta(v2))./(c^2)) + ...
        (1/(2*c^2*gr))*((u1 - u2constant).^2) - E0;
    
    % Giving f one parameter Z = [u1,v1,v2,E0]
    f_oneparam =@(V) f(u2constant, V(1), V(2), E0);
    
    % Optimizer to solve for (v1,v2) values given seed
    f_options = optimset('Display','off','Algorithm','levenberg-marquardt');
    v_star = fsolve(f_oneparam, v_seed, f_options);

    % Update list of (a,v1,v2)
    v_list(counter,1:3) = [a,v_star];

    % Compute Lyapunov exponents of branch?
    if Lyap_bool == true
        % --- Initial conditions for vortices ---
        % Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
        w1_0 = (u2constant) + 1i*(v_star(1)); % positive vortex
        w2_0 = (u2constant) + 1i*(v_star(2)); % negative vortex
        
        % Real and imaginary parts of isothermal coords
        u1_0 = real(w1_0);
        v1_0 = imag(w1_0);
        
        u2_0 = real(w2_0);
        v2_0 = imag(w2_0);

        % Compute MLE
        tf = 1000; % LONGER TIME ---> more accurate MLE
        display = 0; % MLE computation plots off
        params = {N,q,p,c,r,a,R,cap,phi,theta,gr,Dginv};
        y0 = [u1_0, u2_0, v1_0 ,v2_0];
        MLE = Fast_MLE(delta, ExpFac, y0, tf, display, params);

        % Update Lyapunov list of (a, MLE)
        MLE_list(counter,1:2) = [a, MLE];
        fprintf('Counter: %.f and MLE: %f\n',counter,MLE);
    end

    % Update counter
    counter = counter + 1;
end

%%% Make the MLE vs eccentricity plot
%
%figure (1)
%plot(sqrt(1 - (R./MLE_list(:,1)).^2), MLE_list(:,2),'-','LineWidth',2)
%grid on
%xlabel('Eccentricity, $\varepsilon$','Interpreter','latex','FontSize', fs)
%ylabel('MLE, $\lambda_1$','Interpreter','latex','FontSize', fs)
%ax = gca;
%ax.FontSize = fs - 1;


%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
