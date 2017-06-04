function est=est_error(simul_lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,eps1)
ni=size(xi,2);
nb=size(xb,2);
%Matriz de colocación y vector de términos independientes. 
L=matriz(rbf,Lrbf,xi,xb,ctr,eps1);
for i=1:ni
    f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
    f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
%Calculamos los coeficientes de la solución. 
coef=L\f;
%Rejilla de puntos para evaluar el error.
%Generación de puntos de evaluación
neval=40;
v=linspace(0,1,neval);
eval_bd=[[v ;zeros(1,neval)],[v; ones(1,neval)],[zeros(1,neval-2);v(2:neval-1)],[ones(1,neval-2);v(2:neval-1)]];
l1=size(eval_bd,2);
[X,Y]=meshgrid(v(2:end-1),v(2:end-1));
eval_int=[X(:)';Y(:)'];
l2=size(eval_int,2);

for i=1:l1
    aprox(i)=feval('aprox',rbf,ctr,coef,eps1,eval_bd(:,i));
    fun(i)=feval(simul_frontera,eval_bd(1,i),eval_bd(2,i)); %Valor de la solucin real en cada punto.
    error(i)=abs(aprox(i)-fun(i))./max([1 abs(fun(i))]); %Error en el punto. 
end
for i=1:l2
    aprox(i)=feval('aprox',Lrbf,ctr,coef,eps1,eval_int(:,i));
    fun(i)=feval(simul_lu,eval_int(1,i),eval_int(2,i)); %Valor de la solucin real en cada punto.
    error(l1+i)=abs(aprox(i)-fun(i))./max([1 abs(fun(i))]); %Error en el punto. 
end
est=norm(error)/sqrt(l1+l2);
