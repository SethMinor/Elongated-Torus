function [] = plotwrapped(phi, theta, thelinewidth, xlims, ylims, tol, colore)
% Plot data wrapped on the surface of a torus nicely.

% phi = horizontal (toroidal) data
% theta = vertical (poloidal) data 
% (this functions expects that length(phi)=length(theta))

% thelinewidth = LineWidth parameter for plotting (float > 0)
% colore = line segment color (e.g. 'k')
% xlims, ylims = window for plots (e.g. [-2 2])

% tol = tolerance (subtracted from 2pi) for plotting consecutive points;
%       tol = 0 means that any consecutive points separated by less than
%       2pi will be connected by a straight line segment

% create wrapped data
phimod = wrapToPi(phi);
thetamod = wrapToPi(theta);

M = length(phimod);

hold on
for m = 1:M
    if (m == 1) % plot starting position
        if (abs(phimod(m)-phimod(m+1)) < (2*pi-tol)) &&...
                (abs(thetamod(m)-thetamod(m+1)) < (2*pi-tol))
            % in this case, the plot condition is satisfied
            plot(phimod(m:m+1),thetamod(m:m+1),...
                'LineWidth',thelinewidth,'Color',colore)
        end
        plot(phimod(m),thetamod(m),'ok','MarkerSize',8)
    elseif (m > 1) && (m < M) % do regular plotting
        if (abs(phimod(m)-phimod(m+1)) < (2*pi-tol)) &&...
                (abs(thetamod(m)-thetamod(m+1)) < (2*pi-tol))
            % in this case, the plot condition is satisfied
            plot(phimod(m:m+1),thetamod(m:m+1),...,
                'LineWidth',thelinewidth,'Color',colore)
        end
    else % plot finishing position
        plot(phimod(m),thetamod(m),...
            'ok','MarkerFaceColor','k','MarkerSize',8)
    end
end
hold off
grid on
xlim(xlims)
ylim(ylims)
end
