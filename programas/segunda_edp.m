%Ejemplo para la segunda EDP
clear all
format short e
%Tomamos los puntos de colocación interiores iguales a los centros.
%Centros de la frontera fuera del dominio. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=289; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
dyrbf='dyimq' %Simulador de la derivada de la RBF respecto de y. 
simul_u='u2'; %Función simulador de la solución.
simul_lu='lu2'; %Operador de laplace de la función simulador. 
simul_frontera='u2_cond_frontera'; %Simulador de las condiciones de frontera. 
epsilon=1; %Parámetro de forma de las RBF. 
%Generamos los centros.
v=linspace(0,1,sqrt(N));
cont=1;
for i=1:length(v)
    for j=1:length(v)
        ctr(1,cont)=v(i);
        ctr(2,cont)=v(j);
        cont=cont+1;
    end
end

%Puntos interiores y sacamos los centros de la frontera del dominio. 
cont1=1;
cont2=1;
cont3=1;
for k=1:N
    if ctr(1,k)==0
        xb2(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(1,k)=-0.1;
    elseif ctr(1,k)==1
        xb2(:,cont1)=ctr(:,k);
        cont1=cont1+1;
        ctr(1,k)=1.1;
    elseif ctr(2,k)==0
        xb1(:,cont2)=ctr(:,k);
        cont2=cont2+1;
        ctr(2,k)=-0.1;
    elseif ctr(2,k)==1
        xb1(:,cont2)=ctr(:,k);
        cont2=cont2+1;
        ctr(2,k)=1.1;
        
    else
        xi(:,cont3)=ctr(:,k); %puntos de colocación en el interior.
        cont3=cont3+1;
    end  
end
figure(1)
plot(xi(1,:),xi(2,:),'*r')
hold on
plot(xb1(1,:),xb1(2,:),'*b')
hold on
plot(xb2(1,:),xb2(2,:),'*b')

hold off

figure(2)
plot(ctr(1,:),ctr(2,:),'*r')
title('Distribución de los centros')

% Generamos la matriz de colocación 
ni=size(xi,2); %Número de puntos interiores.
nb1=size(xb1,2);
nb2=size(xb2,2);
%Número de puntos en la frontera. 
%Calculamos las matrices de distancias entre los puntos del interior y los
%centros y entre los puntos de la frontera y los centros. 
for j=1:N
    for i=1:ni
        di(i,j)=norm(xi(:,i)-ctr(:,j));
    end
    for i=1:nb1
        db1(i,j)=norm(xb1(:,i)-ctr(:,j));
        dy(i,j)=xb1(2,i)-ctr(2,j);
    end
    for i=1:nb2
        db2(i,j)=norm(xb2(:,i)-ctr(:,j));
         
    end
    
    
end
%Matriz de colocación
for i=1:ni
    for j=1:N
        L(i,j)=feval(Lrbf,di(i,j),epsilon); 
    end
end
for i=1:nb1
    for j=1:N
        L(ni+i,j)=feval(dyrbf,db1(i,j),dy(i,j),epsilon);
    end
end
for i=1:nb2
    for j=1:N
        L(ni+nb1+i,j)=feval(rbf,db2(i,j),epsilon);
    end
end
        

%Vector de términos independientes. 
for i=1:ni
f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb1
f(ni+i,1)=feval(simul_frontera,xb1(1,i),xb1(2,i));
end
for i=1:nb2
f(ni+nb1+i,1)=feval(simul_frontera,xb2(1,i),xb2(2,i));
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
        rbf_eval(j)=feval(rbf,d(j),epsilon);
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