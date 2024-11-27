clear; clc;

d = 2; 

rng(100);
R = diag([1,4]);
U = orth(rand(d,d));
A = U' * R * U ;
A = A ./ vecnorm(A);

c = 1;
fun = @(x) c*d-c*sum(cos(3*pi*x),2)+diag(x*A*x');

K  = 400; lambda = 1/sqrt(d)/2; 

% three fixed initial iterates, as shown in Fig. !
x1three =[0 1.414; 1 -1; 1 1];

% three randomly selected initial iterates
% x1three = 2*rand(3,d)-1;
% x1three = sqrt(d)*x1three./sqrt(sum(x1three.^2,2));

rho1 = 0.93; rho2 = 0.95; rho3 = 0.97;
n = 5;

XTrace1 = app(fun,x1three(1,:),K,lambda,rho1,n);
XTrace2 = app(fun,x1three(2,:),K,lambda,rho2,n);
XTrace3 = app(fun,x1three(3,:),K,lambda,rho3,n);
XTrace01=[x1three(1,:);XTrace1];
XTrace02=[x1three(2,:);XTrace2];
XTrace03=[x1three(3,:);XTrace3];

figure(1)
set(gcf,'Position',[10/0.277 45/0.277 80/0.277 60/0.277]); % 8X6cm
figure_FontSize=10;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

x = linspace(-1.5,1.5);
y = linspace(-1.5,1.5);
[X,Y] = meshgrid(x,y);
u = length(x);
Z = zeros(u,u);
for i=1:u
    for j=1:u
        z = [X(i,j),Y(i,j)];
        Z(i,j) = fun(z);
    end
end
mesh(X,Y,Z);

xticks([-1.5 0 1.5])
xticklabels({'-1.5','0','1.5'})
yticks([-1.5 0 1.5])
yticklabels({'-1.5','0','1.5'})

figure(2)
set(gcf,'Position',[100/0.277 45/0.277 80/0.277 60/0.277]); % 8X6cm
set(gca,'Position',[.15 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

contour(X,Y,Z,10);
xticks([-1.5 0 1.5])
xticklabels({'-1.5','0','1.5'})
yticks([-1.5 0 1.5])
yticklabels({'-1.5','0','1.5'})
hold on
rectangle('position',[0-sqrt(2),0-sqrt(2),2*sqrt(2),2*sqrt(2)],...
    'curvature',[1,1],'EdgeColor','r');
figbN = 5;
plot(XTrace01(1:figbN,1),XTrace01(1:figbN,2),'k-o')
plot(XTrace02(1:figbN,1),XTrace02(1:figbN,2),'b--s')
plot(XTrace03(1:figbN,1),XTrace03(1:figbN,2),'m-.d')
hold off

figure(3)
set(gcf,'Position',[190/0.277 45/0.277 80/0.277 60/0.277]);
set(gca,'Position',[.16 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(log10(sum(XTrace1.^2,2)),'k-')
hold on
plot(log10(sum(XTrace2.^2,2)),'b--')
plot(log10(sum(XTrace3.^2,2)),'m-.')
hold off
ylim([-12 0])
yticks([-12 -8 -4 0])

title(sprintf('d=%d',d))
L1 = sprintf('\x03C1=%03.2f, n=%d',rho1,n);
L2 = sprintf('\x03C1=%03.2f, n=%d',rho2,n);
L3 = sprintf('\x03C1=%03.2f, n=%d',rho3,n);
lgnd=legend(L1,L2,L3);
po=get(lgnd,'Position');
set(lgnd,'Position',[po(1)+0.02, po(2)+0.02, po(3), po(4)]); %8X6cm
xlabel('iteration (k)')
ylabel('$$\log_{10}\|x_k-x_*\|_2^2$$','Interpreter','latex');
