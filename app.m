function [XTrace,YTrace,fc] = app(fun,x1,K,lambda,rho,n)
%APP for global optimization
%By Xiaopeng Luo

d  = length(x1);
xk = x1;
XTrace = zeros(K,d); 
YTrace = zeros(K,1);
alpha  = lambda;
beta   = 1;

p  = haltonset(d,'Skip',1);
p  = scramble(p,'RR2');
fc = inf;
for i=1:K
    % Generate n random vector from halton sequence
    p.Skip = mod((i-1)*n+1,2^52-n);
    x  = net(p,n);
    t  = [xk; norminv(x,xk,1/alpha)];

    % Compute function value sequenc
    f  = fun(t); 
    fk = f(1);
    fc = min(fc,min(f));
    f  = f - min(f);
    if std(f)==0
        break;
    end

    % Use averaged asymptotic formula
    wk = 1e-6 / beta^2 + 1 / ( sqrt(mean(f.^2)) + 1/beta^2 );
    f  = f * wk;
    xk = sum(t.*exp(-f),1) ./ sum(exp(-f));

    % Update and record
    alpha = alpha / rho;
    beta  = beta / rho;
    XTrace(i,:) = xk;
    YTrace(i) = fk;
    fprintf('Iter %d - Objective: %d  dis: %f;\n',i,fk,norm(xk));
end

end

