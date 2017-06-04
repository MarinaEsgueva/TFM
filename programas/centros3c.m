function [ECM,N,ctr,xi,xb,d,coef,iter,condic,eps1]=centros3c(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,tol,dmin,sol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elección de la distribución de centros de las RBF.
%Válida en dominios [0,1]x[0,1]
%EN ENTRADA: 
%   simul_Lu: simulador de la función f(x,y) de la EDP con Lu(x,y)=f(x,y)
%   simul_frontera simulador de la función g(x,y) de las condiciones de
%   frontera: u(x,y)=g(x,y)
%   rbf: base de  rbf elegida
%   Lrbf: operador diferencial aplicado a la RBF.
%   xi:centros iniciales en el interior
%   b1:centros iniciales en [0,1]x 0
%   b2:centros iniciales en [0,1]x 1
%   b3:centros iniciales en 0x[0,1]
%   b4: centros iniciales en 1*[0,1]
%   dmin: mínima distancia permitida entre nodos. 
%   tol: tolerancia en las diferencias entre aproximantes
%   sol: solución real.


% EN SALIDA:
%   ECM: ECM cometido en la resolución de la edp
%   N: número de centros utilizado en la aproximación
%   ctr: conjunto de centros
%   xi: conjunto de puntos de colocación en el interior
%   xb:conjunto de puntos de colocación en la frontera
%   d: menor separación entre centros. 
%   coef: coeficientes de la aproximación
%   iter: número de veces que se ha refinado la rejilla.
%   condic: condicionamiento de la matriz.
%   eps1: parámetro de forma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Incialización de valores. 
ni=size(xi,2); %número de puntos en el interior.
xb=unique([b1 b2 b3 b4]','rows')'; %Se eliminan los posibles puntos repetidos en las esquinas. 
nb=size(xb,2); %Número de puntos en la frontera.
N=ni+nb; %Número de centros
ctr=[xi xb];
%Parámetros para el cálculo del parámetro de forma óptimo. 
minim=0;
maxim=10;
paso=0.2;

d=min([diff(sort(unique(ctr(1,:)))) diff(sort(unique(ctr(2,:)))) ]); %Cota inferior de la distancia ente puntos. 
[eps1,~]=eps_optimo(simul_Lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,minim,maxim,paso); %parámetro de forma óptimo para el conjunto inicial de centros. 

cte=eps1*d; %constante para controlar el condicionamiento. 
ind=true; %Inicializamos para entrar en el bucle while. 
iter=1;
%Se entra en el buble while si hay algún punto test en el que la diferencia
%de aproximaciones es mayor que la toleracia y los centros están lo
%suficientemente separados. 
while (any (ind)) && d>dmin
    %Desplazamos las puntos iniciales interiores
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
      
    %Matriz de distancias para x e y
    for j=1:N
        for i=1:ni
            di(i,j)=norm(xi(:,i)-ctr(:,j));
            di2(i,j)=norm(yi(:,i)-ctr2(:,j));
            
        end
        for i=1:nb
            db(i,j)=norm(xb(:,i)-ctr(:,j));
            db2(i,j)=norm(yb(:,i)-ctr2(:,j));
        end

    end
%Matrices de colocación. 
    for i=1:ni
        Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1); 
        Mi(i,1:N)=feval(Lrbf,di2(i,1:N),eps1);
        
    end
    for i=1:nb
        Lb(i,1:N)=feval(rbf,db(i,1:N),eps1); 
        Mb(i,1:N)=feval(rbf,db2(i,1:N),eps1);
    end
    %Vector de términos independientes. 
    for i=1:ni
        fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
        gi(i,1)=feval(simul_Lu,yi(1,i),yi(2,i));
    end
    for i=1:nb
        fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
        gb(i,1)=feval(simul_frontera,yb(1,i),yb(2,i));
    end
    %Calculamos los coeficientes de la solución. 
    coef=[Li;Lb]\[fi;fb];
    coef2=[Mi;Mb]\[gi;gb];
    condic=cond([Li;Lb]);
    
    %Rejilla de puntos en los que se evalua la diferencia entre soluciones
    
    %Rejilla refinada de puntos interiores. 
    z_x=[0 coord_x]+0.5*dx;
    z_y=[0 coord_y]+0.5*dy;
    [X,Y]=meshgrid(z_x,z_y);
    zi=[X(:)';Y(:)'];
    
    z1=[b1(1,1:end-1)+0.5*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    z2=[b2(1,2:end)-0.5*diff(b2(1,:));ones(1,size(b2,2)-1)];
    z3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.5*diff(b3(2,:))];
    z4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.5*diff(b4(2,:))];
    
    
    %Evaluamos la diferencia entre soluciones  en la rejilla z. 
    ind_int=zeros(1,size(zi,2));
    for i=1:size(zi,2)
        aprox1_int(i)=feval('aprox',rbf,ctr,coef,eps1,zi(:,i));
        aprox2_int(i)=feval('aprox',rbf,ctr2,coef2,eps1,zi(:,i));
    end
    error_int=abs(aprox1_int-aprox2_int)./(1+aprox1_int);
    ind_int(error_int>tol)=1;
    
    %Actualización del conjunto de puntos del interior. 
    xi=[xi zi(:,ind_int==1)];
    ni=size(xi,2);
    
    
    ind_b1=zeros(1,size(z1,2));
    for i=1:size(z1,2)
        aprox1_b1(i)=feval('aprox',rbf,ctr,coef,eps1,z1(:,i));
        aprox2_b1(i)=feval('aprox',rbf,ctr2,coef2,eps1,z1(:,i));
    end
    error_b1=abs(aprox1_b1-aprox2_b1)./(1+aprox1_b1);
    ind_b1(error_b1>tol)=1;
    b1=[b1 z1(:,ind_b1==1)];
    %Reordenar los puntos:
    b1=[sort(b1(1,:));zeros(1,size(b1,2))];
    
    ind_b2=zeros(1,size(z2,2));
    for i=1:size(z2,2)
        aprox1_b2(i)=feval('aprox',rbf,ctr,coef,eps1,z2(:,i));
        aprox2_b2(i)=feval('aprox',rbf,ctr2,coef2,eps1,z2(:,i));
    end
    error_b2=abs(aprox1_b2-aprox2_b2)./(1+aprox1_b2);
    ind_b2(error_b2>tol)=1;
    b2=[b2 z2(:,ind_b2==1)];
    b2=[sort(b2(1,:));ones(1,size(b2,2))];
    
    ind_b3=zeros(1,size(z3,2));
    for i=1:size(z3,2)
        aprox1_b3(i)=feval('aprox',rbf,ctr,coef,eps1,z3(:,i));
        aprox2_b3(i)=feval('aprox',rbf,ctr2,coef2,eps1,z3(:,i));
    end
    error_b3=abs(aprox1_b3-aprox2_b3)./(1+aprox1_b3);
    ind_b3(error_b3>tol)=1;
    b3=[b3 z3(:,ind_b3==1)];
    b3=[zeros(1,size(b3,2));sort(b3(2,:))];
   
    ind_b4=zeros(1,size(z4,2));
    for i=1:size(z4,2)
        aprox1_b4(i)=feval('aprox',rbf,ctr,coef,eps1,z4(:,i));
        aprox2_b4(i)=feval('aprox',rbf,ctr2,coef2,eps1,z4(:,i));
    end
    error_b4=abs(aprox1_b4-aprox2_b4)./(1+aprox1_b4);
    ind_b4(error_b4>tol)=1;
    b4=[b4 z4(:,ind_b4==1)];
    b4=[ones(1,size(b4,2));sort(b4(2,:))];
    %Actualizamos puntos de la frontera
    xb=unique([b1 b2 b3 b4]','rows')';
    nb=size(xb,2);
    %Actualizamos centros.
    ctr=[xi xb];
    N=ni+nb;
    ind=[ind_int ind_b1 ind_b2 ind_b3 ind_b4];
    d=d/2; %Cota inferior de la distancia ente puntos.
    %eps1=cte/d; 
    iter=iter+1;
    
end
%Recalcular la matriz de colocación

for j=1:N
    for i=1:ni
        di(i,j)=norm(xi(:,i)-ctr(:,j));    
    end
    for i=1:nb
        db(i,j)=norm(xb(:,i)-ctr(:,j));
    end
end
for i=1:ni
    Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1); 
end
for i=1:nb
    Lb(i,1:N)=feval(rbf,db(i,1:N),eps1);     
end
for i=1:ni
    fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
end
for i=1:nb
    fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
coef=[Li;Lb]\[fi;fb];
condic=cond([Li;Lb]);
%Rejilla fija para calcular el ECM. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];

for i=1:size(ev_points,2)
    aproxim(i)=aprox(rbf,ctr,coef,eps1,ev_points(:,i));
    fun(i)=feval(sol,ev_points(1,i),ev_points(2,i));
    error(i)=abs(aproxim(i)-fun(i))/max([1,abs(fun(i))]);
    diferencias(i)=abs(aproxim(i)-aprox(rbf,ctr2,coef2,eps1,ev_points(:,i)));
end
Z=reshape(error,40,40);
ZZ=reshape(aproxim,40,40);
ZZZ=reshape(diferencias,40,40);
figure(1)
surf(X,Y,Z);
figure(2)
surf(X,Y,ZZ);
figure(3)
plot(ctr(1,:),ctr(2,:),'*')
figure(4)
surf(X,Y,ZZZ);
figure(5)
contour(X,Y,Z)
figure(6)
contour(X,Y,ZZZ)
ECM=1/neval*norm(error);