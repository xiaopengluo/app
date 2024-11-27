function [XTrace,YTrace,best] = de(fun,d,itermax,para)
% Differential Evolution 
% from book 'Nature-Inspired Optimization Algorithms' by Xinshe Yang
% 2022-06-020

XTrace = zeros(itermax,d); 
YTrace = zeros(itermax,1);

n = para(1);    % Population >=4, typically 10 to 25
F = para(2);    % DE parameter - scaling (0.5 to 0.9)
Cr = para(3);   % DE parameter - crossover probability

% Simple bounds [-1,1]^d
Lb = -ones(d,1); 
Ub = ones(d,1);

% Dimension of search variables
d = length(Lb);

% Initialize the population/solutions
for i=1:n
    Sol(i,:) = Lb+(Ub-Lb).*rand(size(Lb));
    Fitness(i) = fun(Sol(i,:));
end
% Find the current best
[fmin,I] = min(Fitness);
best = Sol(I,:); 

% Start the iterations by differential evolution
for iter=1:itermax
    % Obtain donor vectors by permutation
    k1 = randperm(n);   k2 = randperm(n);
    k1sol = Sol(k1,:);  k2sol = Sol(k2,:);
    % Random crossover index/matrix
    K = rand(n,d)<Cr;
    % DE/RAND/1 scheme
    % 正如注释所说，这里采用了 DE-rand-1 的更新策略，还有其他更新策略：DE-best-1 DE-rand/to/best-1
    V = Sol+F*(k1sol-k2sol);
    V = Sol.*(1-K)+V.*K;
    
    % Evaluate new solutions
    for i=1:n
        Fnew = fun(V(i,:));
        % If the solution improves
        if Fnew <= Fitness(i)
            Sol(i,:) = V(i,:);
            Fitness(i) = Fnew;
        end
        % Update the current best
        if Fnew <= fmin
            best = V(i,:); 
            fmin = Fnew;
        end
    end
    XTrace(iter,:) = best;
end