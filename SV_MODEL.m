%%sv model
Lam=0.0233; 
lambda=2.5;
Gam=7.4; 
gamma=4.3;
t1=[0 1 2 3 4 5 6]; 
p_cluster=Lam*exp(-Lam*t1); % ideal exponential pdf
plot(t1, p_cluster,'r*') 
