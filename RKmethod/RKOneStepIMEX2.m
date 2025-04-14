function pt = RKOneStepIMEX2(odeEX,t0,dt, p0,G,L,BT,BCRhs)

Fs{1}=dt*odeEX(t0,p0);

K{1}=p0+BT.aE(2,1)*Fs{1};

K{1}=L\BCRhs(K{1},p0);
Fs{2}=dt*odeEX(t0+BT.cE(2)*dt,K{1});
Gs{1}=dt*G*K{1};
K{2}=p0+BT.aE(3,1)*Fs{1}+BT.aE(3,2)*Fs{2}+ BT.aI(3,2)*Gs{1};
pt=L\BCRhs(K{2},K{1});%pt=K{2};



end