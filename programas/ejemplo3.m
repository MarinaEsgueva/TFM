clear all
sol='u3';
simul_Lu='Lu3';
simul_frontera='u3';
rbf='imq';
Lrbf='Limq';
xi=[0.25 0.25 0.5 0.75 0.75  ;0.25 0.75  0.5  0.25 0.75 ];
b1=[0 0.5 1; 0 0 0];
b2=[0 0.5 1; 1 1 1];
b3=[0 0 0 ;0 0.5 1];
b4=[1 1 1 ;0 0.5 1];
dmin=0.05;
Nmax=50;
indic=1;
[ECM,N,ctr,xi,xb,coef,condic,eps_int,eps_front]=centros5b(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,dmin,sol,Nmax,indic);
 
