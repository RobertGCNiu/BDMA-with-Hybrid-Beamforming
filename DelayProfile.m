 figure
 h=stem(rays(:,3)*10^9, (20*log10(abs(rays(:,1)))+30), '.');
 set(h,'BaseValue',-180);
 xlabel('Delay (ns)')
 ylabel('Ray power (dBm)')
 grid on