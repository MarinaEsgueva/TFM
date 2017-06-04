clear all
format short e
% C�lculo del par�metro de forma �ptimo. 
%Centros de la forntera fuera del dominio.
%En el interior los puntos de colocaci�n y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=169; %N�mero de puntos de colocaci�n. 
neval=40; %Tama�o de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tama�o neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_frontera='u'; %Funci�n simulador.
simul_lu='Lu'; %Operador de laplace de la funci�n simulador. 
%Generamos los centros en una rejilla uniforme. 
v=linspace(0,1,sqrt(N));
cont=1;
for i=1:length(v)
    for j=1:length(v)
        ctr(1,cont)=v(i);
        ctr(2,cont)=v(j);
        cont=cont+1;
    end
end

%Puntos de colocaci�n interiores y de frontera. 
%Sacamos los centros de la frontrera del dominio. 
cont1=1;
cont2=1;
for k=1:N
    if ctr(1,k)==0
        xb(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(1,k)=-0.1;
    elseif ctr(1,k)==1
        xb(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(1,k)=1.1;
    elseif ctr(2,k)==0
        xb(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(2,k)=-0.1;
    elseif ctr(2,k)==1
        xb(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(2,k)=1.1;
        
    else
        xi(:,cont2)=ctr(:,k); %puntos de colocaci�n en el interior.
        cont2=cont2+1;
    end  
end



%Variamos el par�metro epsilon 
cont=1;
for epsilon=0.1:0.2:10
    ep(cont)=epsilon;
    error(cont)=generalizacion_Rippa(simul_lu,simul_frontera,rbf,Lrbf,epsilon,ctr,xi,xb);
    cont=cont+1;
end

semilogy(ep,error)
[estimacion_ecm,n]=min(error);
optimo=ep(n);
