
function [H,a_TX,a_RX]=NearAoD(Num_users,TX_ant_w,RX_ant_w,Num_paths,K)

%Lam=0.0233;

M_k = Num_users/K;
TX_ant_h = TX_ant_w;
RX_ant_h = RX_ant_w;
H=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);  % One user channel
a_TX=zeros(TX_ant_w*TX_ant_h,Num_users,Num_paths); % TX steering vector
a_TX_left = zeros(TX_ant_w*TX_ant_h,Num_users,Num_paths);
a_RX=zeros(RX_ant_w*RX_ant_h,Num_users,Num_paths); % RX steering vector
a_RX_left=zeros(RX_ant_w*RX_ant_h,Num_users,Num_paths); % RX steering vector
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);
% Constructing the channels
for k=1:1:K
    if k==1
        AoD_el_k=pi*rand(1,Num_paths)-pi/2;
        AoD_az_k=2*pi*rand(1,Num_paths);
        AoA_el_k=pi*rand(1,Num_paths)-pi/2;
        AoA_az_k=2*pi*rand(1,Num_paths);
    elseif k==2
        AoD_el_k=AoD_el_k+(k-1)*pi/2;
        AoD_az_k=AoD_az_k+(k-1)*pi/2;
        AoA_el_k=AoA_el_k+(k-1)*pi/2;
        AoA_az_k=AoA_az_k+(k-1)*pi/2;
    end
    
    for u = 1+(k-1)*M_k:1:k*M_k
        AoD_el(u,:)=AoD_el_k+ rand(1,Num_paths)/5;
        AoD_az(u,:)= AoD_az_k +   rand(1,Num_paths)/5;
        AoA_el(u,:)=AoA_el_k + rand(1,Num_paths)/5;
        AoA_az(u,:)=AoA_az_k + rand(1,Num_paths)/5;
        alpha(u,:)=  sqrt(1/Num_paths)*sqrt(1/2)*(randn(1,Num_paths)+1j*randn(1,Num_paths));
        %     alpha(u,:)= Lam*exp(-Lam*(0:20:(Num_paths-1)*20)); % ideal exponential pdf
        %     alpha(u,:) = alpha(u,:) ./sum(alpha);
        
        Temp_Channel=zeros(RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);
        for l=1:1:Num_paths
            a_TX(:,u,l)=transpose(sqrt(1/(TX_ant_w*TX_ant_h))*exp(1j*pi*(ind_TX_w*sin(AoD_az(u,l))*sin(AoD_el(u,l))+ind_TX_h*cos(AoD_el(u,l))) ));
            a_RX(:,u,l)=transpose(sqrt(1/(RX_ant_w*RX_ant_h))*exp(1j*pi*(ind_RX_w*sin(AoA_az(u,l))*sin(AoA_el(u,l))+ind_RX_h*cos(AoA_el(u,l))) ));
            Temp_Channel=Temp_Channel+sqrt((TX_ant_w*TX_ant_h)*(RX_ant_w*RX_ant_h))*alpha(u,l)*a_RX(:,u,l)*a_TX(:,u,l)';
        end
        H(u,:,:)=Temp_Channel;
    end
end

end
