function est=generalizacion_Rippa(simul_lu,simul_frontera,rbf,Lrbf,epsilon,ctr,xi,xb)
ni=size(xi,2); 
nb=size(xb,2);
N=ni+nb;
%Calculamos los coeficientes de la matriz de colocación. 
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
f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
for k=1:ni
    indices=[1:N];
    indices(k)=[]; %Eliminamos el punto k-esimo 
    L_temp=L(indices,indices); %Eliminamos la fila y la columna k-esimas 
    f_temp=f(indices); %Eliminamos la componente k-esima
    coef_temp=L_temp\f_temp; %Calculamos los coeficientes de la aproximación. 
    %Calculamos la Ek
    cont=1;
    for j=indices  
            d=norm(xi(:,k)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            rbf_eval(cont)=feval(Lrbf,d,epsilon);
            cont=cont+1;  
    end
    aprox=rbf_eval*coef_temp;
    error(k)=aprox-feval(simul_lu,xi(1,k),xi(2,k));
end
for k=1:nb
    indices=[1:N];
    indices(k)=[]; %Eliminamos el punto k-esimo 
    L_temp=L(indices,indices); %Eliminamos la fila y la columna k-esimas 
    f_temp=f(indices); %Eliminamos la componente k-esima
    coef_temp=L_temp\f_temp; %Calculamos los coeficientes de la aproximación. 
    %Calculamos la Ek
    cont=1;
    for j=indices  
            d=norm(xb(:,k)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            rbf_eval(cont)=feval(rbf,d,epsilon);
            cont=cont+1;  
    end
    aprox=rbf_eval*coef_temp;
    error(k+ni)=aprox-feval(simul_frontera,xb(1,k),xb(2,k));
end
est=norm(error);
