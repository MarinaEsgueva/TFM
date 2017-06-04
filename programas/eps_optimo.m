function [eps,est]=eps_optimo(simul_lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,minim,max,paso)
N=size(ctr,2);
ni=size(xi,2);
nb=size(xb,2);
%Matriz de distancias
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

%Generación de puntos de evaluación
neval=40;
v=linspace(0,1,neval);
eval_bd=[[v ;zeros(1,neval)],[v; ones(1,neval)],[zeros(1,neval-2);v(2:neval-1)],[ones(1,neval-2);v(2:neval-1)]];
l1=size(eval_bd,2);
[X,Y]=meshgrid(v(2:end-1),v(2:end-1));
eval_int=[X(:)';Y(:)'];
l2=size(eval_int,2);


cont=1;
for epsilon=minim:paso:max
    param(cont)=epsilon;
    %Matriz de colocación.
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
    for i=1:l1
        for j=1:N
            d(j)=norm(eval_bd(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            rbf_eval(j)=feval(rbf,d(j),epsilon);
        end
        aprox(i,1)=rbf_eval*coef; %Calculamos la aproximación con RBF en cada punto
        fun(i,1)=feval(simul_frontera,eval_bd(1,i),eval_bd(2,i)); %Valor de la solucin real en cada punto.
        error(i,1)=abs(aprox(i,1)-fun(i,1)); %Error en el punto. 
    end
    for i=1:l2
        for j=1:N
            d(j)=norm(eval_int(:,i)-ctr(:,j)); %Distancias de los puntos de evaluacion a los centros.
            Lrbf_eval(j)=feval(Lrbf,d(j),epsilon);
        end
        aprox(i,1)=Lrbf_eval*coef; %Calculamos la aproximación con RBF en cada punto
        fun(i,1)=feval(simul_lu,eval_int(1,i),eval_int(2,i)); %Valor de la solucin real en cada punto.
        error(l1+i,1)=abs(aprox(i,1)-fun(i,1)); %Error en el punto. 
    end
    ecm(cont)=norm(error)/(l1+l2); %Error cuadrático medio
    cont=cont+1;
    
end
[est n]=min(ecm);
eps=param(n); 



