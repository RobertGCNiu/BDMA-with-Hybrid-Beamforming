theta = pi/3; theta2 = pi/4; theta3 = pi/5;
n = 0:100; m = 0:60; 
a_T = exp(1j*pi*n*sin(theta));
a_R = exp(1j*pi*m*sin(theta));
a_T2 = exp(1j*pi*n*sin(theta2));
a_R2 = exp(1j*pi*m*sin(theta2));
a_T3 = exp(1j*pi*n*sin(theta3));
a_R3 = exp(1j*pi*m*sin(theta3));

1/6000*a_R*(a_R'*a_T)*a_T2'
1/6000*a_R*(a_R3'*a_T3)*a_T2'