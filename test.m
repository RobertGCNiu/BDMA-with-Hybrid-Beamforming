clear all;
Num_users=16;
N=12;
[H,a_TX,a_RX]=generate_channels(user,N,N,N,N,1);
h1=reshape(H(1,:,:),[N^2,N^2]);
% h2=reshape(H(2,:,:),[N^2,N^2]);

tx1=a_TX(:,1);
rx1=a_RX(:,1);

Frf=a_TX;
Wrf=a_RX;

I=0;
for i=2:user
    I=I+abs(a_RX(:,i)'*h1*a_TX(:,i))^2;
end
SIR= abs(rx1'*h1*tx1)^2/I;
rate=log2(1+SIR)

