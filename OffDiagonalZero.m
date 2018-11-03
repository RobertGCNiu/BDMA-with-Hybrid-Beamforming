function F_odz = OffDiagonalZero(N_RF, N_users,G)
F_1 = pinv(G(1:N_RF,1:N_RF));
F_2 = pinv(G(N_RF+1:N_users, N_RF+1:N_users));
F_odz = blkdiag(F_1, F_2);
end

