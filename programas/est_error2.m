function est=est_error2(simul_lu,simul_frontera,xi,xb,rbf,Lrbf,eps_int,eps_front)
ni=size(xi,2);
nb=size(xb,2);
N=ni+nb;
%Matriz de colocación y vector de términos independientes. 
L1=matriz(rbf,Lrbf,xi,xb,xi,eps_int);
L2=matriz(rbf,Lrbf,xi,xb,xb,eps_front);
for i=1:ni
    f(i,1)=feval(simul_lu,xi(1,i),xi(2,i));
end
for i=1:nb
    f(ni+i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
end
%Calculamos los coeficientes de la solución. 
coef=[L1 L2]\f;
%Rejilla de puntos para evaluar el error.
%Generación de puntos de evaluación
neval=40;
v=linspace(0,1,neval);
eval_bd=[[v;zeros(1,neval)],[v; ones(1,neval)],[zeros(1,neval-2);v(2:neval-1)],[ones(1,neval-2);v(2:neval-1)]];
l1=size(eval_bd,2);
[X,Y]=meshgrid(v(2:end-1),v(2:end-1));
eval_int=[X(:)';Y(:)'];
l2=size(eval_int,2);

for i=1:l1
    aprox(i)=feval('aprox2',rbf,xi,xb,coef,eps_int,eps_front,eval_bd(:,2));
    fun(i)=feval(simul_frontera,eval_bd(1,i),eval_bd(2,i)); %Valor de la solución real en cada punto.
    error(i)=abs(aprox(i)-fun(i))./max([1 abs(fun(i))]); %Error en el punto. 
end
for i=1:l2
    aprox(i)=feval('aprox2',Lrbf,xi,xb,coef,eps_int,eps_front,eval_int(:,i));
    fun(i)=feval(simul_lu,eval_int(1,i),eval_int(2,i)); %Valor de la solucin real en cada punto.
    error(l1+i)=abs(aprox(i)-fun(i))./max([1 abs(fun(i))]); %Error en el punto. 
end
est=norm(error)/sqrt(l1+l2);
