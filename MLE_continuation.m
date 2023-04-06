%% Branching for MLE plot
clear, clc;%, clf;

% Font size, for plotting
fs = 16;

% Initial parameters
a = 13;
a0 = a;
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
da = 0.05;
counter = 1;

for delta_a = 0:da:1
    % New semi-major axis a
    a = a0 - delta_a;

    % Initial seed
    if counter == 1
        V = v_seed;
    else
        V = v_star;
    end
    
    % --- Solving for isothermal coordinates ---
    % Boundary condition
    phi_0 = 0;
    
    % Numerical integration using ode45
    theta_span = [0, pi];
    options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);
    [theta_raw,phi_raw] = ode15s(@(theta,phi) odefcn(theta,phi,a,R,r), theta_span, phi_0, options);
    
    % Give these hosses odd symmetry 
    theta_raw = [-flip(theta_raw(2:end)); theta_raw];
    phi_raw = [-flip(phi_raw(2:end)); phi_raw];
    
    % Quick plot to verify isothermal solutions look good
    v0 =@(theta) 2*r*atan(sqrt((R-r)/(R+r))*tan(theta/2));
    
    % Define (phi,theta) to (u,v) map
    vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp'); % phi = f(theta)
    u =@(phi) c*phi;
    v =@(theta) -c*vhelper(theta);
    
    % Define gl and gr
    gl = min(imag(phi_raw));
    gr = max(imag(phi_raw));
    
    % Defining the (u,v) to (phi,theta) map
    ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
    phi =@(u) u/c;
    theta =@(v) ginv(-UVwrap(v, [-c*gr, c*gr])/c);
    
    % --- Other parameters ---
    % Gamma quantity
    gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
    
    % Local scale factor
    lambda =@(phi,theta) gamma(phi,theta)./c;
    
    % Jacobi theta function parameters
    cap = 20;     % truncation error
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

    % Update counter
    counter = counter + 1;
end


%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
