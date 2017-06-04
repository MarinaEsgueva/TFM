function [ECM,N,ctr,xi,xb,coef,condic,eps_int,eps_front]=centros4b(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,dmin,sol,Nmax,indic)
%indic: indicador de parámetro de forma fijo o variable 
% 0: fijo
% 1: distinguir interior y frontera
%Incialización de valores. 
ni=size(xi,2); %número de puntos en el interior.
xb=unique([b1 b2 b3 b4]','rows')'; %Se eliminan los posibles puntos repetidos en las esquinas. 
nb=size(xb,2); %Número de puntos en la frontera.
N=ni+nb; %Número de centros
ctr=[xi xb];

%Parámetros para el cálculo del parámetro de forma óptimo. 
minim=0.1;
maxim=10;
paso=0.2;


%Generamos la rejilla de posibles nuevos centros:
v=[0:dmin:1];
[X,Y]=meshgrid(v,v);
new_ctr=[X(:)';Y(:)'];
M=size(new_ctr,2);
plot(new_ctr(1,:),new_ctr(2,:),'*')
%Eliminamos de la rejilla los centros
ind=zeros(1,M);
for i=1:N
    for j=1:M
        if norm(new_ctr(:,j)-ctr(:,i))<0.05
            ind(j)=1;
        end
    end
end
 new_ctr(:,ind==1)=[];
 M=size(new_ctr,2);
 
 
%Rejilla fija para calcular el ECM. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];

    
while N<Nmax
    
    %Parámetro de forma óptimo 
    [eps_int,fval,exitflag]=fminbnd(@(eps1)est_error(simul_Lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,eps1),minim,maxim);
    if indic==0
        eps_front=eps_int; 
    else
        [eps_front,fval,exitflag]=fminbnd(@(eps_front)est_error2(simul_Lu,simul_frontera,xi,xb,rbf,Lrbf,eps_int,eps_front),minim,eps_int);
    end
    
    %Desplazamiento de los centros y puntos interiores. 
    coord_x=sort(unique(xi(1,:))); 
    coord_y=sort(unique(xi(2,:)));
    dx=diff([0 coord_x 1]);
    dy=diff([0 coord_y 1]);
    d1=min([dx dy]);
    yi=xi+0.3*d1*ones(2,ni);
    
    %Puntos de la frontera. Hay que actualizar lado a lado 
    bb1=[b1(1,1:end-1)+0.3*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    bb2=[b2(1,2:end)-0.3*diff(b2(1,:));ones(1,size(b2,2)-1)];
    bb3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.3*diff(b3(2,:))];
    bb4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.3*diff(b4(2,:))];
    yb=[bb1 bb2 bb3 bb4];
    ctr2=[yi yb];
    
    L1=matriz(rbf,Lrbf,xi,xb,xi,eps_int);
    L2=matriz(rbf,Lrbf,xi,xb,xb,eps_front);
    M1=matriz(rbf,Lrbf,yi,yb,yi,eps_int);
    M2=matriz(rbf,Lrbf,yi,yb,yb,eps_front);
    %Vector de términos independientes. 
    for i=1:ni
        f(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
        g(i,1)=feval(simul_Lu,yi(1,i),yi(2,i));
    end
    for i=1:nb
        f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
        g(ni+i,1)=feval(simul_frontera,yb(1,i),yb(2,i));
    end
    %Calculamos los coeficientes de la solución. 
    coef=[L1 L2]\f;
    coef_old=coef;
    coef2=[M1 M2]\g;
    
    
    %Diferencia entre aproximaciones en los posibles nuevos centros
    diferencia=[];
    for i=1:M
        aprox1=feval('aprox2',rbf,xi,xb,coef,eps_int,eps_front,new_ctr(:,i));
        aprox2=feval('aprox2',rbf,yi,yb,coef2,eps_int,eps_front,new_ctr(:,i));
        diferencia(i)=abs(aprox1-aprox2)/max([1,abs(aprox1)]);
    end
    [m,n]=max(diferencia);
    nuevo=new_ctr(:,n);
    new_ctr(:,n)=[];
    M=M-1;
    xi_old=xi; %Guardamos el último conjunto de centros para representar las diferencias. 
    xb_old=xb; 
    % Comprobamos a que conjunto de centros iniciales pertenece
    if abs(nuevo(2))<1.e-8
        b1=[b1 nuevo];
        b1=[sort(b1(1,:));zeros(1,size(b1,2))];
        nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(2)-1)<1.e-8
        b2=[b2 nuevo];
        b2=[sort(b2(1,:));ones(1,size(b2,2))];
        nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(1))<1.e-8
        b3=[b3 nuevo];
        b3=[zeros(1,size(b3,2));sort(b3(2,:))];
        nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(1)-1)<1.e-8
        b4=[b4 nuevo];
        b4=[ones(1,size(b4,2));sort(b4(2,:))];
        nb=nb+1;
        xb=[xb nuevo];
    else
        xi=[xi nuevo];
        ni=ni+1;
    end
    N=N+1;
    
    ctr=[ctr nuevo]; % Añadir el nuevo centro
    
    
end
%Actualización de la matriz de colocación 
L1=matriz(rbf,Lrbf,xi,xb,xi,eps_int);
L2=matriz(rbf,Lrbf,xi,xb,xb,eps_front);
%Vector de términos independientes. 
for i=1:ni
    f(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
end
for i=1:nb
    f(i+ni,1)=feval(simul_frontera,xb(1,i),xb(2,i));    
end
%Calculamos los coeficientes de la aproximación. 
coef=[L1 L2]\f;
condic=cond([L1 L2]);


for i=1:size(ev_points,2)
    aproxim(i)=feval('aprox2',rbf,xi,xb,coef,eps_int,eps_front,ev_points(:,i));
    fun(i)=feval(sol,ev_points(1,i),ev_points(2,i));
    err(i)=abs(aproxim(i)-fun(i))/max([1,abs(fun(i))]);
    dif(i)=abs(feval('aprox2',rbf,xi_old,xb_old,coef_old,eps_int,eps_front,ev_points(:,i))-feval('aprox2',rbf,yi,yb,coef2,eps_int,eps_front,ev_points(:,i)))/max([1,abs(feval('aprox2',rbf,xi_old,xb_old,coef_old,eps_int,eps_front,ev_points(:,i)))]);
end
Z=reshape(err,40,40);
ZZ=reshape(aproxim,40,40);
ZZZ=reshape(dif,40,40);
figure(1)
surf(X,Y,Z);
title('Error')
figure(2)
surf(X,Y,ZZ);
figure(3)
plot(ctr(1,:),ctr(2,:),'*',nuevo(1),nuevo(2),'r*')
hold on 
contour (X,Y,ZZZ)
figure(4)
surf(X,Y,ZZZ);
title('Diferencias')
ECM=1/neval*norm(err);