clear;clc
theta1_t = pi/3;
theta2_t = pi/4;
theta1_r = pi/5;
theta2_r = pi/6;
N_R = 64;
n_R = 0:63;
times = 0;
interference_r_match = 1/N_R^2*abs( (sum(exp(1j.*n_R.*pi.*(sin(theta1_r)-sin(theta1_r)))))).^2;
interference_r_no_match = 1/N_R^2*abs( (sum(exp(1j.*n_R.*pi.*(sin(theta1_r)-sin(theta2_r)))))).^2;
N_T_range = 12:12:144;
L = 2;
for N_T = N_T_range
    times = times+1;
n = 0:N_T-1;

interference_v(times) = 1/N_T^2*abs( (sum(exp(1j.*n.*pi.*(sin(theta1_t)-sin(theta2_t)))))).^2;
end

inter1 =  interference_v;
inter2 =  interference_r_no_match * interference_v;
plot(N_T_range, inter1);
hold on
plot(N_T_range, inter2);
plot(N_T_range, ones(length(N_T_range),1))