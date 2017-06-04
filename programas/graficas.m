%Representar las soluciones aproximadas gráficamentes
neval=50;
grid=linspace(0,1,neval)
[X,Y]=meshgrid(grid,grid);
ev_points=[X(:)',Y(:)'];
rbf='imq'
coef=coef;
centros=ctr;
N=size(centros,2);
eps1=1;
for i=1:neval^2
    for j=1:N
        rbf_eval(j)=feval(rbf,norm(ev_points(:,i)-centros(:,j)),eps1);
    end
    aprox(i)=rbf_eval*coef;
end
aprox=reshape(aprox,neval,neval);
surf(X,Y,aprox)