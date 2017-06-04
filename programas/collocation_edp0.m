clear all
format short e
%Caso EDP homogenea. 
%Tomamos los puntos de colocación interiores iguales a los centros.
%Centros de la forntera fuera del dominio
N=289; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_u='u'; %Función simulador.
simul_part='u_part2'; %Simulador de la solución particular.
epsilon=5; %Parámetro de forma de las RBF. 
ind=2;
[ctr,xi,xb,ev_points]=generacion_puntos(N,neval,ind);
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
    for j=1:N
        L(i,j)=feval(Lrbf,di(i,j),epsilon); 
    end
end
for i=1:nb
    for j=1:N
        L(ni+i,j)=feval(rbf,db(i,j),epsilon);
    end
end
%Vector de términos independientes. 
for i=1:ni
f(i,1)=0; 
end
for i=1:nb
f(ni+i,1)=feval(simul_u,xb(1,i),xb(2,i))-feval(simul_part,xb(1,i),xb(2,i));
end
%Calculamos los coeficientes de la solución. 
coef=L\f;

%Evaluación de la función en los puntos de evaluacion. 
for i=1:neval^2
    for j=1:N
        d(j)=norm(ev_points(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
        rbf_eval(j)=feval(rbf,d(j),epsilon);
    end
    aprox(i,1)=rbf_eval*coef +feval(simul_part,ev_points(1,i),ev_points(2,i)); %Calculamos la aproximación con RBF en cada punto
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