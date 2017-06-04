function est=est_ecm(simul_lu,simul_frontera,rbf,Lrbf,epsilon,xi,xb,ctr,N)
ni=size(xi,2); %Número de puntos interiores
nb=size(xb,2); %Número de puntos en la frontera
%Resolvemos la EDP con todos los puntos:
%Generación de la matriz de distancias en el interior y en la frontera
for j=1:N
    for i=1:ni %Filas correspondientes al interior
        di(i,j)=norm(xi(:,i)-ctr(:,j));
    end
    for i=1:nb %Filas correspondientes a la frontera
        db(i,j)=norm(xb(:,i)-ctr(:,j));
    end
    
end
for i=1:ni %Términos intependientes del interior
    f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb %Terminos independientes de la frontera. 
    f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
%Matriz de colocación
for j=1:N
        for i=1:ni
            L(i,j)=feval(Lrbf,di(i,j),epsilon); 
        end
   
        for i=1:nb
            L(ni+i,j)=feval(rbf,db(i,j),epsilon);
        end
end
%Cálculo de los coeficientes de la aproximación. 
coef=L\f;
%Generación de los puntos de evaluación. 
neval=40;
grid=linspace(0,1,neval);
[xe,ye]=meshgrid(grid);
ev_points=[xe(:) ye(:)]';
%Evaluación de la función en los puntos de evaluacion. 
for i=1:neval^2
    for j=1:N
        d=norm(ev_points(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
        rbf_eval(j)=feval(rbf,d,epsilon);
    end
    aprox(i,1)=rbf_eval*coef; %Calculamos la aproximación con RBF en cada punto 
end
%Resolvemos un problema eliminando un punto del conjunto de centros y de
%puntos interiores. 
for k=1:N %Eliminamos un centro cada vez
    %Vector de indices que identifica con que centros y puntos de
    %colocación trabajamos en cada iteración. 
    indices=[1:N];
    indices(k)=[]; %Eliminamos el punto k-esimo 
    L_temp=L(indices,indices); %Eliminamos la fila y la columna k-esimas 
    f_temp=f(indices); %Eliminamos la componente k-esima
    coef_temp=L_temp\f_temp; %Calculamos los coeficientes de la aproximación. 
    %Calculamos la nueva aproximación en cada uno de los puntos de
    %colocación. 
    for i=1:neval^2
        cont=1;
        rbf_eval=[];
        for j=indices  
            d=norm(ev_points(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            rbf_eval(cont)=feval(rbf,d,epsilon);
            cont=cont+1;  
        end
        aprox_2(i,1)=rbf_eval*coef_temp; %Calculamos la nueva aproximación con RBF en cada punto 
    end
 error(k)=(1/neval)*norm(aprox-aprox_2);       
end
est=norm(error) 
    



    