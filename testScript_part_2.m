Pw=Pparts.*wsave;
Zw=Zparts.*wsave;
Nw=Nparts.*wsave;
Dw=Dparts.*wsave;
gw=gammaparts.*wsave;

P = sum(Pw,2);
Z = sum(Zw,2);
N = sum(Nw,2);
D = sum(Dw,2);
G = sum(gw,2);

figure
hold on
plot(1:548,P(:,1),'-r')
plot(1:548,Z(:,1),'-b')
plot(1:548,N(:,1),'-m')
plot(1:548,D(:,1),'-g')
plot(1:548,G(:,1),'-k')
title('marine ecosystem');
legend('phytoplankton','zooplankton','nutrients','detritus','\gamma');
xlabel('time')
ylabel('{grams of carbon}/{meter^3}');
hold off