
clear all
format short e
%Primera EDP 
%Centros de la forntera fuera del dominio.
%En el interior los puntos de colocación y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
%Segunda estrategia de variación
N=1089; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_u='u'; %Función simulador.
simul_lu='Lu'; %Operador de laplace de la función simulador. 
epsilon_min=10; %Parámetro de forma de las RBF. 
epsilon_max=1;
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

%Puntos de colocación interiores y de frontera. 
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
        xi(:,cont2)=ctr(:,k); %puntos de colocación en el interior.
        cont2=cont2+1;
    end  
end
figure(1)
plot(xi(1,:),xi(2,:),'*r')
hold on
plot(xb(1,:),xb(2,:),'*b')
hold off

figure(2)
plot(ctr(1,:),ctr(2,:),'*r')
title('Distribución de los centros')

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
        L(i,j)=feval(Lrbf,di(i,j),sqrt(epsilon_min^2*(epsilon_max^2/epsilon_min^2)^((j-1)/(N-1)))); 
    end
end
for i=1:nb
    for j=1:N
        L(ni+i,j)=feval(rbf,db(i,j),sqrt(epsilon_min^2*(epsilon_max^2/epsilon_min^2)^((j-1)/(N-1))));
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
    for j=1:N
        d(j)=norm(ev_points(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
        rbf_eval(j)=feval(rbf,d(j),sqrt(epsilon_min^2*(epsilon_max^2/epsilon_min^2)^((j-1)/(N-1))));
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