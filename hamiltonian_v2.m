function [energy,classic,curve,quantum] = hamiltonian_v2(U,V,N,q,p,c,r,a,R,cap,phi,theta,gr)
% This function computes the total energy of a system of vortices.

% INPUTS
% W = U + iV = wortex positions vector (in isothermal coorindates)
% N = number of vortices
% q = vector of vortex charges (should sum to 0)
% r,R,c = torus parameters
% p = nome of first Jacobi theta function (should satisfy |p| < 1)
% cap = truncation term of sum in first Jacobi theta function

% Define ye olde complex coordinate
W = U + 1i*V;

% Gamma
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);

% Local scale factor
L =@(u,v) gamma(phi(u),theta(v))/c;

% Initialize energy contributions
classic = 0; curve = 0; quantum = 0;

% Loop over each vortex pair
for n = 1:N
    % Curvature contribution
    curve = curve + log(L(U(:,n),V(:,n)));

    for m = 1:N
        % quantum contribution
        quantum = quantum + (1/(2*c^2*gr)).*q(n).*q(m).*U(:,n).*U(:,m);

        if (m ~= n)
            % classical contribution
            thetus = (W(:,n) - W(:,m))./(2*c);
            classic = classic -...
                (q(n).*q(m).*log(abs(jacobitheta1(thetus,p,cap))));
        end
    end
end
energy = classic + curve + quantum;
end