%Script para elección de centros.
clear all
sol='u';
simul_Lu='Lu';
simul_frontera='u';
rbf='imq';
Lrbf='Limq';
xi=[0.25 0.5 0.5 0.75; 0.75 0.1 0.9 0.75];
b1=[0 0.25 0.5 0.75 1; zeros(1,5)];
b2=[0 0.5 1; ones(1,3)];
b3=[0 0 0;0 0.5 1];
b4=[1 1 1;0 0.5 1];
dmin=0.05;
Nmax=50;
indic=1;
[ECM,N,ctr,xi,xb,coef,cond,eps_int,eps_front]=centros5b(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,dmin,sol,Nmax,indic);

