clear all
format short e
% Cálculo del parámetro de forma óptimo. 
%Centros de la forntera fuera del dominio.
%En el interior los puntos de colocación y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=289; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='gaussiana'; 
Lrbf='Lgaussiana'; %Operador de laplace de la RBF. 
simul_frontera='u3'; %Función simulador.
simul_u='u3';
simul_lu='Lu3'; %Operador de laplace de la función simulador. 
ind=2;
[ctr,xi,xb,ev_points,xe,ye]=generacion_puntos(N,neval,ind);

% Generamos la matriz de colocación 
ni=size(xi,2); %Número de puntos interiores.
nb=size(xb,2);  %Número de puntos en la frontera. 
%Calculamos las matrices de distancias entre los puntos del interior y los
%centros y entre los puntos de la frontera y los centros. 
for j=1:N
    for i=1:ni
        di(i,j)=norm(xi(:,i)-ctr(:,j));
    end
    for i=1:nb
        db(i,j)=norm(xb(:,i)-ctr(:,j));
    end
    
end
for i=1:ni
f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
%Variamos el parámetro epsilon 
cont=1;
for epsilon=0.1:0.2:10
    ep(cont)=epsilon;
    for i=1:ni
        for j=1:N
            L(i,j)=feval(Lrbf,di(i,j),epsilon); 
        end
    end
    for i=1:nb
        for j=1:N
            L(ni+i,j)=feval(rbf,db(i,j),epsilon);
        end
    end
 
    %Calculamos los coeficientes de la solución. 
    coef=L\f;
    %Evaluación de la función en los puntos de evaluacion. 
    for i=1:neval^2
        for j=1:N
            d(j)=norm(ev_points(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            rbf_eval(j)=feval(rbf,d(j),epsilon);
        end
        aprox(i,1)=rbf_eval*coef; %Calculamos la aproximación con RBF en cada punto
        fun(i,1)=feval(simul_u,ev_points(1,i),ev_points(2,i)); %Valor de la solucin real en cada punto.
        error(i,1)=abs(aprox(i,1)-fun(i,1)); %Error en el punto. 
    end
    ecm(cont)=norm(error)/neval; %Error cuadrático medio
    est_Rippa(cont)=generalizacion_Rippa(simul_lu,simul_frontera,rbf,Lrbf,epsilon,ctr,xi,xb);
    cont=cont+1;
end

semilogy(ep,[ecm;est_Rippa])
[error,n]=min(ecm);
[est,m]=min(est_Rippa);
optimo=ep(n);
optimo_Rippa=ep(m);
error_est=ecm(m)