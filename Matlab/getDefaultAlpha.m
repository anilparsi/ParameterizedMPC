function alpha = getDefaultAlpha(Ac,Bc,Q,R)
% Estimate the decay rate of Laguerre basis functions using the closed loop
% LQR solution
[~,~,E] = lqr(Ac,Bc,Q,R);
alpha = max(real(E));
end