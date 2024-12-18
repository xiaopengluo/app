function [XTrace,YTrace,fc] = app(fun,x1,K,lambda,rho,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d  = length(x1);
xk = x1;
XTrace = zeros(K,d); 
YTrace = zeros(K,1);
alpha  = lambda;

p  = haltonset(d);
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

    % Use averaged asymptotic formula
    if mean(f)~=0
        f = f / mean(f);
    end
    xk = sum(t.*exp(-f),1) ./ sum(exp(-f));

    % Update and record
    alpha = alpha / rho;
    XTrace(i,:) = xk;
    YTrace(i) = fk;
    fprintf('Iter %d - Objective: %d;\n',i,fk);
end
end

