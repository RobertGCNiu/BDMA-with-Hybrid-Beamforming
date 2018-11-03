function [y,tor_optimal] = fixPointIter(q,func,heff, P_max, noise_k, tol,MaxIter)
% Version: 1.0 written by jbb0523 @2016-08-21
    if nargin < 5
        MaxIter = 100;
    end
	if nargin < 4
        tol = 1e-3;
    end
    xn =  q;
   % fprintf('Iter  0: %16.14f\n',q);
    xnp1 =  func(q,1, heff,noise_k,P_max);
  %  fprintf('Iter  1: %16.14f\n',xnp1);
    criterion = abs(xnp1-xn);
    xn = xnp1;
    Iter = 1;
    while(criterion>tol)
        for u = 1: length(q)
            q_k(u) = func(q,u, heff,noise_k,P_max);
        end
        criterion = sum(abs(q_k-q));
        q = q_k;
        %fprintf('Iter  1: %16.14f\n',criterion);
        if Iter>=MaxIter
            break;
        end
    end
    y = q;
    [~, tor_optimal] = func(y,1, heff,noise_k,P_max);
end
