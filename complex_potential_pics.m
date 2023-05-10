% Pictures of the complex potential F
% Stream function Re(F), phase contours Im(F)
clear, clc;

% Fontsize, for plotting
fs = 14;

R = 12;
r = 9;
c = sqrt(R^2 - r^2);
p = exp(-pi*r*c);
cap = 20;

%jacobitheta1(w,p,cap)
w1 = -1 + 0*1i;
w2 = +1 - 0.5*1i;

% Complex potential
F =@(w) log(jacobitheta1((w-w1)./(2*c),p,cap) ./ jacobitheta1((w-w2)./(2*c),p,cap))...
    - w.*real(w1 - w2)./(2*pi*r*c);
% Stream function
figure (1)
ax1 = subplot(2,2,1);
[X,Y] = meshgrid(-2*pi:0.1:2*pi, -2*pi:0.1:2*pi);
Chi = real(F(X+1i*Y));
s1 = surf(X,Y,Chi);
colormap(ax1, parula);
camlight('left')
shading flat;
%s1.EdgeColor = 'none';

%xlabel('','Interpreter','latex','FontSize',fs)
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
title('Stream Function, $\chi$','Interpreter','latex','FontSize',fs)

ax2 = subplot(2,2,2);
Phi = imag(F(X+1i*Y));
s2 = surf(X,Y,Phi);
colormap(ax2,hot);
camlight('left')
shading flat;

xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
title('Phase Field, $\Phi$','Interpreter','latex','FontSize',fs)


ax3 = subplot(2,2,3:4);
surf(X,Y,Chi-4);
hold on
countour_levels = [-3, -2.5:0.4:2.5, 3];
contour(X,Y,Phi,countour_levels,'-w','LineWidth',0.5)
contour(X,Y,Chi,countour_levels,'-k','LineWidth',1)
hold off
colormap(ax3, parula);
clim([-8, 0])
%camlight('right')
shading flat;
view([0,0,90])

xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
%legend('','$\chi$','$\Phi$','Interpreter','latex','FontSize',fs,'Color',[211,211,211]./256)
