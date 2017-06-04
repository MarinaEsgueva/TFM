clear all
format short e
%Primera EDP
%Parámetro de forma diferente en la frontera y en el interior.
%En el interior se hace variar de forma aleatoria alrededor de la media del
%interior y analogamente en la frontera. 
%Centros de la frontera fuera del dominio.
%En el interior los puntos de colocación y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=1089; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
ind=2;
rbf='imq'; 

Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_u='u3'; %Función simulador.
simul_lu='Lu3'; %Operador de laplace de la función simulador. 
epsilon_interior=2.1; %Parámetro de forma en el interior  
epsilon_frontera=1.3; %Parámetro de forma en la frontera
[ctr_bd,ctr_ind,xi,xb,ev_points,xe,ye]=generacion_puntos2(N,neval,ind);
% Generamos la matriz de colocación 
ni=size(xi,2); %Número de puntos y centros interiores.
nb=size(xb,2);  %Número de puntos y centros en la frontera. 
%Generamos los vectores de parámetros de forma en el interior y en la
%frontera
eps_interior=epsilon_interior*(0.5*ones(ni,1)+rand(ni,1));
eps_frontera=epsilon_frontera*(0.5*ones(nb,1)+rand(nb,1));




%Calculamos las matrices de distancias entre los puntos del interior y los
%centros y entre los puntos de la frontera y los centros. 

%Distancias a los centros del interior
for j=1:ni
    for i=1:ni
        dii(i,j)=norm(xi(:,i)-ctr_int(:,j));
    end
    for i=1:nb
        dbi(i,j)=norm(xb(:,i)-ctr_int(:,j));
    end
    
end
for j=1:nb
    for i=1:ni
        dib(i,j)=norm(xi(:,i)-ctr_bd(:,j));
    end
    for i=1:nb
        dbb(i,j)=norm(xb(:,i)-ctr_bd(:,j));
    end
    
end
for i=1:ni
    for j=1:ni
        L(i,j)=feval(Lrbf,dii(i,j),eps_interior(j)); 
    end
    for j=1:nb
        L(i,ni+j)=feval(Lrbf,dib(i,j),eps_frontera(j)); 
    end
end
for i=1:nb
    for j=1:ni
        L(ni+i,j)=feval(rbf,dbi(i,j),eps_interior(j));
    end
    for j=1:nb
        L(ni+i,ni+j)=feval(rbf,dbb(i,j),eps_frontera(j));
    end
end
%Vector de términos independientes. 
for i=1:ni
f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
f(ni+i,1)=feval(simul_u,xb(1,i),xb(2,i));
end
%Calculamos los coeficientes de la solución. 
coef=L\f;
%Generamos la rejilla de los puntos de evaluacion. 
grid=linspace(0,1,neval);
[xe,ye]=meshgrid(grid);
ev_points=[xe(:) ye(:)]';
%Evaluación de la función en los puntos de evaluacion. 
for i=1:neval^2
    for j=1:ni
        d(j)=norm(ev_points(:,i)-ctr_int(:,j)); %Distancias de los puntos de evaluacion a los centros interiores.
        rbf_eval(j)=feval(rbf,d(j),eps_interior(j));
    end
    for j=1:nb
        d(ni+j)=norm(ev_points(:,i)-ctr_bd(:,j)); %Distancias de los puntos de evaluacion a los centros interiores.
        rbf_eval(ni+j)=feval(rbf,d(ni+j),eps_frontera(j));
    end
    aprox(i,1)=rbf_eval*coef; %Calculamos la aproximación con RBF en cada punto
    fun(i,1)=feval(simul_u,ev_points(1,i),ev_points(2,i)); %Valor de la solucin real en cada punto.
    error(i,1)=abs(aprox(i,1)-fun(i,1)); %Error en el punto. 
end
ecm=norm(error)/neval %Error cuadrático medio
cond(L) %Condicionamiento de la matriz de colocación. 
aprox=reshape(aprox,neval,neval);
fun=reshape(fun,neval,neval);
error=reshape(error,neval,neval);
figure(3)
surf(xe,ye,fun)
title('Solución real')
figure(4)
surf(xe,ye,aprox)
title('Solución aproximada')
figure(5)
surf(xe,ye,error)
title('Distribución del error')