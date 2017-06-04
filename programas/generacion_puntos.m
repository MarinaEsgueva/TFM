function[ctr,xi,xb,ev_points,xe,ye]=generacion_puntos(N,neval,ind)
%Genera los centros, los puntos de colocaci�n y los puntos donde se evalua
%el error.
%Input:
%N: n�mero de puntos de centros y de puntos de colocaci�n.
%neval: tama�o de la rejilla de puntos de evaluaci�n (es de tama�o
%neval*neval)
% ind= toma los siguientes valores
%      1: centros en la frontera.
%      2: centros de la frontera fuera del dominio.
%Output: 
%ctr: centros de las RBF
%xi: puntos de colocaci�n interiores
%xb: puntos de colocaci�n en la frontera
%ev_points: puntos en los que se evalua el error.
%xe e ye: coordenadas de la rejilla. 
%Generamos los centros en una rejilla uniforme. 
v=linspace(0,1,sqrt(N));
nv=length(v);
%Puntos de colocaci�n del borde.
xb=[[zeros(1,nv);v],[ones(1,nv);v],[v(2:end-1);zeros(1,nv-2)],[v(2:end-1);ones(1,nv-2)]];
if ind==2
   ctr=[[-0.1*ones(1,nv);v],[1.1*ones(1,nv);v],[v(2:end-1);-0.1*ones(1,nv-2)],[v(2:end-1);1.1*ones(1,nv-2)]];   
else
    ctr=xb;
end
v(1)=[];
v(end)=[];
[X,Y]=meshgrid(v,v);
ctr=[ctr,[X(:)'; Y(:)']];
xi=[X(:)';Y(:)'];

%Generamos la rejilla de los puntos de evaluacion. 
grid=linspace(0,1,neval);
[xe,ye]=meshgrid(grid);
ev_points=[xe(:) ye(:)]';
    
figure(1)
plot(xi(1,:),xi(2,:),'*r')
hold on
plot(xb(1,:),xb(2,:),'*b')
hold off

figure(2)
plot(ctr(1,:),ctr(2,:),'*r')
title('Distribuci�n de los centros')