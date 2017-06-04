function [ECM,N,ctr,xi,xb,d,coef,iter,condic,eps1]=centros3(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,tol,sol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elección de la distribución de centros de las RBF.
%Comienza con N0 nodos distribuidos en una rejilla y va añadiendo centros
%en puntos donde el error sea mayor que una tolerancia dada. 
%Válida en dominios [0,1]x[0,1]
%EN ENTRADA: 
%   simul_Lu: simulador de la función f(x,y) de la EDP con Lu(x,y)=f(x,y)
%   simul_frontera simulador de la función g(x,y) de las condiciones de
%   frontera: u(x,y)=g(x,y)
%   rbf: base de  rbf elegida
%   Lrbf: operador diferencial aplicado a la RBF.
%   xi:centros iniciales en el interior
%   b1:centros iniciales en [0,1]x 0
%   b2:centros iniciales en [0,1]x 1
%   b3:centros iniciales en 0x[0,1]
%   b4: centros iniciales en 1*[0,1]
%   tol: tolerancia en las diferencias entre aproximantes
%   sol: solución real.


% EN SALIDA:
%   ECM :ECM cometido en la resolución de la edp
%   N: número de centros utilizado en la aproximación
%   ctr: conjunto de centros
%   xi: conjunto de puntos de colocación en el interior
%   xb:conjunto de puntos de colocación en la frontera
%   d: menor separación entre centros. 
%   coef: coeficientes de la aproximación
%   iter: número de veces que se ha refinado la rejilla.
%   condic: condicionamiento de la matriz.
%   eps1: parámetro de forma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Incialización de valores. 
ni=size(xi,2); %número de puntos en el interior.
xb=unique([b1 b2 b3 b4]','rows')';
nb=size(xb,2); %Número de puntos en la frontera.
N=ni+nb;
ctr=[xi b1 b2 b3 b4]; 
minim=0;
maxim=10;
paso=0.2;
%d=min([diff(unique(ctr(1,:))) diff(unique(ctr(1,:))) ]); %Cota inferior de la distancia ente puntos. 
[eps1,~]=eps_optimo(simul_Lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,minim,maxim,paso);
%cte=eps1*d;
ind=true; %Inicializamos para entrar en el bucle while. 
iter=1;

while (any (ind)) 
    %Desplazamos las puntos iniciales interiores
    coord_x=unique(xi(1,:));
    coord_y=unique(xi(2,:));
    dx=diff([0 coord_x 1]);
    dy=diff([0 coord_y 1]);
    d=min([dx dy]);
    yi=xi+0.3*d*ones(2,ni);
   
    %Puntos de la frontera. Hay que actualizar lado a lado 
    bb1=[b1(1,1:end-1)+0.3*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    bb2=[b2(1,2:end)-0.3*diff(b2(1,:));ones(1,size(b2,2)-1)];
    bb3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.3*diff(b3(2,:))];
    bb4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.3*diff(b4(2,:))];
    yb=[bb1 bb2 bb3 bb4];
   
    ctr2=[yi yb];
    
    
    %Matriz de distancias para x e y
    for j=1:N
        for i=1:ni
            di(i,j)=norm(xi(:,i)-ctr(:,j));
            di2(i,j)=norm(yi(:,i)-ctr2(:,j));
            
        end
        for i=1:nb
            db(i,j)=norm(xb(:,i)-ctr(:,j));
            db2(i,j)=norm(yb(:,i)-ctr2(:,j));
        end

    end
%Matriz de colocación. 
    for i=1:ni
        Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1); 
        Mi(i,1:N)=feval(Lrbf,di2(i,1:N),eps1);
        
    end
    for i=1:nb
        Lb(i,1:N)=feval(rbf,db(i,1:N),eps1); 
        Mb(i,1:N)=feval(rbf,db2(i,1:N),eps1);
    end
    %Vector de términos independientes. 
    for i=1:ni
        fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
        gi(i,1)=feval(simul_Lu,yi(1,i),yi(2,i));
    end
    for i=1:nb
        fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
        gb(i,1)=feval(simul_frontera,yb(1,i),yb(2,i));
    end
    %Calculamos los coeficientes de la solución. 
    coef=[Li;Lb]\[fi;fb];
    coef2=[Mi;Mb]\[gi;gb];
    condic=cond([Li;Lb]);

    %Rejilla de puntos en los que se evalua la diferencia entre soluciones
    
    %Rejilla refinada de puntos interiores. 
    z_x=[0 coord_x]+0.5*dx;
    z_y=[0 coord_y]+0.5*dy;
    [X,Y]=meshgrid(z_x,z_y);
    zi=[X(:)';Y(:)'];
    
    z1=[b1(1,1:end-1)+0.5*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    z2=[b2(1,2:end)-0.5*diff(b2(1,:));ones(1,size(b2,2)-1)];
    z3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.5*diff(b3(2,:))];
    z4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.5*diff(b4(2,:))];
    
    
    %Evaluamos la diferencia entre soluciones  en la rejilla z. 
    ind_int=zeros(1,size(zi,2));
    for i=1:size(zi,2)
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(zi(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(zi(:,i)-ctr2(:,j)),eps1);
        end
        if abs(rbf_eval1*coef-rbf_eval2*coef2)/max([1 abs(rbf_eval1*coef)])>tol
            ind_int(i)=1; %Añadimos el punto al conjunto de puntos de colocación. 
        end
    end
    %Actualización del conjunto de puntos del interior. 
    xi=[xi zi(:,ind_int==1)];
    ni=size(xi,2);
    
    
    ind_b1=zeros(1,size(z1,2));
    for i=1:size(z1,2);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z1(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z1(:,i)-ctr2(:,j)),eps1);
        end
        if abs(rbf_eval1*coef-rbf_eval2*coef2)/max([1 abs(rbf_eval1*coef)])>tol
            ind_b1(i)=1;
        end
    end
    b1=[b1 z1(:,ind_b1==1)];
    
    ind_b2=zeros(1,size(z2,2));
    for i=1:size(z2,2);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z2(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z2(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)/max([1 abs(rbf_eval1*coef)])>tol
            ind_b2(i)=1;
        end
    end
    b2=[b2 z2(:,ind_b2==1)];
    
    ind_b3=zeros(1,size(z3,2));
    for i=1:size(z3,2);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z3(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z3(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)/max([1 abs(rbf_eval1*coef)])>tol
            ind_b3(i)=1;
        end
    end
    
    b3=[b3 z3(:,ind_b3==1)];
    
    ind_b4=zeros(1,size(z4,2));
    for i=1:size(z4,2);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z4(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z4(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)/max([1 abs(rbf_eval1*coef)])>tol
            ind_b4(i)=1;
        end
    end
    b4=[b4 z4(:,ind_b4==1)];
    %Actualizamos puntos de la frontera
    xb=unique([b1 b2 b3 b4]','rows')';
    nb=size(xb,2);
    %Actualizamos centros.
    ctr=[xi xb];
    N=ni+nb
    ind=[ind_int ind_b1 ind_b2 ind_b3 ind_b4]; 
    %d=min([diff(unique(ctr(1,:))) diff(unique(ctr(1,:))) ]); %Cota inferior de la distancia ente puntos.
    %eps1=cte/d; 
    iter=iter+1;
    
end

%Rejilla fija para calcular el ECM. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];

for i=1:size(ev_points,2)
    for j=1:N
        rbf_eval(j)=feval(rbf,norm(ev_points(:,i)-ctr(:,j)),eps1);
    end
    aprox(i)=rbf_eval*coef;
    error(i)=abs(feval(sol,ev_points(1,i),ev_points(2,i))-rbf_eval*coef)/max([1,abs(feval(sol,ev_points(1,i),ev_points(2,i)))]);
end
Z=reshape(error,40,40);
ZZ=reshape(aprox,40,40);
figure(1)
surf(X,Y,Z);
figure(2)
surf(X,Y,ZZ);
figure(3)
plot(ctr(1,:),ctr(2,:),'*')
ECM=1/neval*norm(error);