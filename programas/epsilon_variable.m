clear all
format short e
%Primera EDP
%Parámetro de forma diferente en la frontera y en el interior
%Centros de la frontera fuera del dominio.
%En el interior los puntos de colocación y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=81; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_u='u3'; %Función simulador.
simul_lu='Lu3'; %Operador de laplace de la función simulador. 
epsilon_interior=0.3; %Parámetro de forma en el interior  
epsilon_frontera=0.1; %Parámetro de forma en la frontera
const_bd=0; %Constante para variar el parámtro alrededor de la media en la frontera.
const_int=0; % Constante par variar el parámetro de forma alrededor de la media en el interior.
%Generamos los centros del interior y los de la frontera 
v=linspace(0,1,sqrt(N));
cont=1;
for i=1:length(v)
    ctr_bd(1,cont)=-0.1;
    ctr_bd(2,cont)=v(i);
    xb(1,cont)=0;
    xb(2,cont)=v(i);
    cont=cont+1;
end
for i=1:length(v)
    ctr_bd(1,cont)=1.1;
    ctr_bd(2,cont)=v(i);
    xb(1,cont)=1;
    xb(2,cont)=v(i);
    cont=cont+1;
end
for i=2:(length(v)-1)
    ctr_bd(1,cont)=v(i);
    ctr_bd(2,cont)=-0.1;
    xb(1,cont)=v(i);
    xb(2,cont)=0;
    cont=cont+1;
end
for i=2:(length(v)-1)
    ctr_bd(1,cont)=v(i);
    ctr_bd(2,cont)=1.1;
    xb(1,cont)=v(i);
    xb(2,cont)=1;
    cont=cont+1;
end
cont=1;
for i=2:(length(v)-1)
    for j=2:(length(v)-1)
        ctr_int(1,cont)=v(i);
        ctr_int(2,cont)=v(j);
        cont=cont+1;
    end
end
xi=ctr_int; %En el interior los puntos de colocación coinciden con los centros. 


figure(1)
plot(xi(1,:),xi(2,:),'*r')
hold on
plot(xb(1,:),xb(2,:),'*b')
title('Distribución de los puntos de colocación')
hold off

figure(2)
plot(ctr_int(1,:),ctr_int(2,:),'*r')
hold on 
plot(ctr_bd(1,:),ctr_bd(2,:),'*r')
title('Distribución de los centros')

% Generamos la matriz de colocación 
ni=size(xi,2); %Número de puntos y centros interiores.
nb=size(xb,2);  %Número de puntos y centros en la frontera. 
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
        L(i,j)=feval(Lrbf,dii(i,j),epsilon_interior/(sqrt(1+const_int*(-1)^j))); 
    end
    for j=1:nb
        L(i,ni+j)=feval(Lrbf,dib(i,j),epsilon_frontera/(sqrt(1+const_bd*(-1)^j))); 
    end
end
for i=1:nb
    for j=1:ni
        L(ni+i,j)=feval(rbf,dbi(i,j),epsilon_interior/(sqrt(1+const_int*(-1)^j)));
    end
    for j=1:nb
        L(ni+i,ni+j)=feval(rbf,dbb(i,j),epsilon_frontera/(sqrt(1+const_bd*(-1)^j)));
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
        rbf_eval(j)=feval(rbf,d(j),epsilon_interior/(sqrt(1+const_int*(-1)^j)));
    end
    for j=1:nb
        d(ni+j)=norm(ev_points(:,i)-ctr_bd(:,j)); %Distancias de los puntos de evaluacion a los centros interiores.
        rbf_eval(ni+j)=feval(rbf,d(ni+j),epsilon_frontera/(sqrt(1+const_bd*(-1)^j)));
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