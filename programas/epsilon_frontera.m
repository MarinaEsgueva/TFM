clear all
format short e
%Calcula el parámetro óptimo en la frontera fijado el epsilon en el
%interior. 
%Primera EDP 
%Centros de la forntera fuera del dominio.
%En el interior los puntos de colocación y los centros coinciden. 
%Distribuidos en una rejilla en el intervalo [0,1]
N=81; %Número de puntos de colocación. 
neval=40; %Tamaño de la rejilla de los puntos donde se mide el error.
          %La rejilla es de tamaño neval*neval.
rbf='imq'; 
Lrbf='Limq'; %Operador de laplace de la RBF. 
simul_u='u3'; %Función simulador.
simul_lu='Lu3'; %Operador de laplace de la función simulador. 
epsilon_interior=0.3; %Parámetro de forma de las RBF en el interior. 
ind=2;
[ctr_bd,ctr_ind,xi,xb,ev_points,xe,ye]=generacion_puntos2(N,neval,ind);
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
% Terminos independientes
for i=1:ni
f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
f(ni+i,1)=feval(simul_u,xb(1,i),xb(2,i));
end
%Generamos la rejilla de los puntos de evaluacion. 
grid=linspace(0,1,neval);
[xe,ye]=meshgrid(grid);
ev_points=[xe(:) ye(:)]';
%Variamos el parametro de forma en la fontera.
cont=1;
for ep_front=0.1:0.1:epsilon_interior
    ep(cont)=ep_front;
    


    for i=1:ni
        for j=1:ni
            L(i,j)=feval(Lrbf,dii(i,j),epsilon_interior); 
        end
        for j=1:nb
            L(i,ni+j)=feval(Lrbf,dib(i,j),ep_front); 
        end
    end
    for i=1:nb
        for j=1:ni
            L(ni+i,j)=feval(rbf,dbi(i,j),epsilon_interior);
        end
        for j=1:nb
            L(ni+i,ni+j)=feval(rbf,dbb(i,j),ep_front);
        end
    end

    %Calculamos los coeficientes de la solución. 
    coef=L\f;

    %Evaluación de la función en los puntos de evaluacion. 
    for i=1:neval^2
        for j=1:ni
            d(j)=norm(ev_points(:,i)-ctr_int(:,j)); %Distancias de los puntos de evaluacion a los centros interiores.
            rbf_eval(j)=feval(rbf,d(j),epsilon_interior);
        end
        for j=1:nb
            d(ni+j)=norm(ev_points(:,i)-ctr_bd(:,j)); %Distancias de los puntos de evaluacion a los centros interiores.
            rbf_eval(ni+j)=feval(rbf,d(ni+j),ep_front);
        end
        aprox(i,1)=rbf_eval*coef; %Calculamos la aproximación con RBF en cada punto
        fun(i,1)=feval(simul_u,ev_points(1,i),ev_points(2,i)); %Valor de la solucin real en cada punto.
        error(i,1)=abs(aprox(i,1)-fun(i,1)); %Error en el punto. 
    end
    ecm(cont)=norm(error)/neval; %Error cuadrático medio
    cont=cont+1;
end
semilogy(ep,ecm)
[error,n]=min(ecm);
optimo=ep(n);