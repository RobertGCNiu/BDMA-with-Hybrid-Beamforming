function F_odz = OffDiagonalZero(N_RF, N_users,F_full)

F_Blockpart = blkdiag(ones(N_RF),ones(N_users-N_RF));
F_odz=F_full.* F_Blockpart;

end

