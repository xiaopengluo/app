clear; clc;

d = 100; 

a = 2; % mean
b = 4; % std

rng(100);
pd = makedist('Normal');
t = truncate(pd,0,inf);
r = random(t,d,1);
r = a + r*b ;
R = diag(r) ;
rng(100);
U = orth(rand(d,d));
A = U' * R * U ;
A = A ./ vecnorm(A);
c = 1;

fun = @(x) c*d-c*sum(cos(3*pi*x),2)+diag(x*A*x');
K  = 6000; lambda = 1/sqrt(d); 

rng('shuffle');
x1three = 2*rand(2,d)-1;
x1three = sqrt(d)*x1three./sqrt(sum(x1three.^2,2));

rho1 = 0.9940; 
n = 25;

XTrace = app(fun,x1three(1,:),K,lambda,rho1,n);

Kde = 12000; para = [20 0.95 0.04];

XTrace0 = de(fun,d,Kde,para);

% various scale factors (F)
F = [.9 .95 1 1.05];
XTrace1 = de(fun,d,Kde,[20 F(1) para(3)]);
XTrace2 = de(fun,d,Kde,[20 F(2) para(3)]);
XTrace3 = de(fun,d,Kde,[20 F(3) para(3)]);
XTrace4 = de(fun,d,Kde,[20 F(4) para(3)]);

% various crossover parameters (Cr)
Cr = [.02 .03 .04 .05];
XTrace5 = de(fun,d,Kde,[20 para(2) Cr(1)]);
XTrace6 = de(fun,d,Kde,[20 para(2) Cr(2)]);
XTrace7 = de(fun,d,Kde,[20 para(2) Cr(3)]);
XTrace8 = de(fun,d,Kde,[20 para(2) Cr(4)]);

figure(1)
figure_FontSize=10;
set(gcf,'Position',[20/0.277 25/0.277 80/0.277 60/0.277]); % 8X6cm
set(gca,'Position',[.16 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(para(1)*(1:Kde),log10(sum(XTrace0.^2,2)),'k-')
ylim([-12 8])
yticks([-12 -8 -4 0 4 8])
hold on
plot(n*(1:K),log10(sum(XTrace.^2,2)),'b--')
hold off

title(sprintf('d=%d',d))
L1 = sprintf('DE:n=%d,F=%03.2f,Cr=%03.2f',para(1),para(2),para(3));
L2 = sprintf('APP: \x03C1=%05.4f, n=%d',rho1,n);
lgnd=legend(L1,L2);
po=get(lgnd,'Position');
set(lgnd,'Position',[po(1)+0.02, po(2)+0.02, po(3), po(4)]); %8X6cm
xlabel('nfe')
ylabel('$$\log_{10}\|x_k-x_*\|_2^2$$','Interpreter','latex');

figure(2)
figure_FontSize=10;
set(gcf,'Position',[110/0.277 25/0.277 80/0.277 60/0.277]); % 8X6cm
set(gca,'Position',[.16 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(para(1)*(1:Kde),log10(sum(XTrace1.^2,2)),'k:')
% xlim([0 15e4])
ylim([-12 8])
yticks([-12 -8 -4 0 4 8])
hold on
plot(20*(1:Kde),log10(sum(XTrace2.^2,2)),'b--')
plot(20*(1:Kde),log10(sum(XTrace3.^2,2)),'m-.')
plot(20*(1:Kde),log10(sum(XTrace4.^2,2)),'g-')
hold off

title('DE: various scale factors')
L1 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),F(1),para(3));
L2 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),F(2),para(3));
L3 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),F(3),para(3));
L4 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),F(4),para(3));
lgnd=legend(L1,L2,L3,L4);
po=get(lgnd,'Position');
set(lgnd,'Position',[po(1)+0.02, po(2)+0.02, po(3), po(4)]); %8X6cm
xlabel('nfe')
ylabel('$$\log_{10}\|x_k-x_*\|_2^2$$','Interpreter','latex');

figure(3)
figure_FontSize=10;
set(gcf,'Position',[200/0.277 25/0.277 80/0.277 60/0.277]); % 8X6cm
set(gca,'Position',[.16 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(para(1)*(1:Kde),log10(sum(XTrace5.^2,2)),'k:')
% xlim([0 15e4])
ylim([-12 8])
yticks([-12 -8 -4 0 4 8])
hold on
plot(20*(1:Kde),log10(sum(XTrace6.^2,2)),'b--')
plot(20*(1:Kde),log10(sum(XTrace7.^2,2)),'m-.')
plot(20*(1:Kde),log10(sum(XTrace8.^2,2)),'g-')
hold off

title('DE: various crossover parameters')
L1 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),para(2),Cr(1));
L2 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),para(2),Cr(2));
L3 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),para(2),Cr(3));
L4 = sprintf('n=%d,F=%03.2f,Cr=%03.2f',para(1),para(2),Cr(4));
lgnd=legend(L1,L2,L3,L4);
po=get(lgnd,'Position');
set(lgnd,'Position',[po(1)+0.02, po(2)+0.02, po(3), po(4)]); %8X6cm
xlabel('nfe')
ylabel('$$\log_{10}\|x_k-x_*\|_2^2$$','Interpreter','latex');