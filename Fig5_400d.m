clear; clc;

d = 400; 

a = 1; % Mean
b = 2; % Std
rng(100); % Fix the random number seed
pd = makedist('Normal'); % Create normal distribution object
t = truncate(pd,0,inf); % Truncate the normal distribution object
r = random(t,d,1); % Generate a random dx1 vector from t
r = a + r*b ; % Eigenvalue vector
R = diag(r) ; % Eigenvalue matrix
rng(100); % Fix the random number seed
U = orth(rand(d,d)); % Obtain an orthonormal basis
A = U' * R * U ; % Compute A
A = A ./ vecnorm(A); % Column normalization

c = 1;

fun = @(x) c*d-c*sum(cos(3*pi*x),2)+diag(x*A*x');

K  = 20000; lambda = 1/sqrt(d); 

rng('shuffle');
x1three = 2*rand(3,d)-1;
x1three = sqrt(d)*x1three./sqrt(sum(x1three.^2,2));


rho1 = 0.9986; rho2 = 0.9988; rho3 = 0.9990; 

n = 45;

XTrace1 = app(fun,x1three(1,:),K,lambda,rho1,n);
XTrace2 = app(fun,x1three(2,:),K,lambda,rho2,n);
XTrace3 = app(fun,x1three(3,:),K,lambda,rho3,n);

figure(1)
figure_FontSize=10;
set(gcf,'Position',[100/0.277 25/0.277 80/0.277 60/0.277]); % 8X6cm
set(gca,'Position',[.16 .18 .78 .72]); % 8X6cm
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

plot(log10(sum(XTrace1.^2,2)),'k-')
xlim([0 K])
ylim([-12 4])
yticks([-12 -8 -4 0 4])
hold on
plot(log10(sum(XTrace2.^2,2)),'b--')
plot(log10(sum(XTrace3.^2,2)),'m-.')
hold off

title(sprintf('d=%d',d))
L1 = sprintf('\x03C1=%05.4f, n=%d',rho1,n);
L2 = sprintf('\x03C1=%05.4f, n=%d',rho2,n);
L3 = sprintf('\x03C1=%05.4f, n=%d',rho3,n);
lgnd=legend(L1,L2,L3);
po=get(lgnd,'Position');
set(lgnd,'Position',[po(1)+0.02, po(2)+0.02, po(3), po(4)]); %8X6cm
xlabel('iteration (k)')
ylabel('$$\log_{10}\|x_k-x_*\|_2^2$$','Interpreter','latex');