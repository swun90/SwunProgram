function Diamond
pdeSolv(1)
pdeSolv(0)
end

function pdeSolv(varargin)
% initialize
global a1 a2 a3 b1 b2 b3 c1 c2 d
a1 = 0.2;
a2 = 0.7;
a3 = 0.7;
b1 = 1.5;
if(nargin<1)
    b2 = 1; %b2 - mean flow
elseif (nargin == 1)
    b2 = varargin{1};
end
b3 = 1;
c1 = 1;
c2 = 0.5; %c2 - neoclassic
d  = 1;

% dt = 0.001;
% Q = 0.01 .* t;

t = [0, 2/0.01];
% y0 = [0.02; 0.01; 0.001];
y0 = [0.02;0.01;0.001];
[t,y] = ode45(@(t,y) ppme(t,y), t, y0);
Q = t*0.01;
%%
figure;
plot(Q, y(:,1), 'r-','LineWidth',2);
hold on
plot(Q, y(:,2), 'b-.','LineWidth',1.5);
plot(Q, y(:,3)/5, 'g--','LineWidth',1.5);
set(gca,'FontSize',14);
l1=legend('$\varepsilon$','$V_{ZF}$','N/5');
set(l1,'interpreter','latex',...
    'location','best',...
    'FontSize',14);
xlabel('Q');
% ylim([0,1.5]);
title(sprintf('Evolution of pp model (b2=%d)',b2));
print(gcf,'-dpng',sprintf('Evol_pp_model-(b2=%d)',b2));
pause(0.5)
close
end