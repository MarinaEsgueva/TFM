%Representación gráfica de Lu

%Generación de la rejilla.

x=linspace(0,1,40);
[X,Y]=meshgrid(x,x);
points=[X(:)';Y(:)'];
fun='u3';
fun=feval(fun,points(1,:),points(2,:));
Z=reshape(fun,40,40);
surf(X,Y,Z)