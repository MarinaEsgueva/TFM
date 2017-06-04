%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLACION CON RBF, 1-dim, PROCESO ADAPTATIVO
% Ref: Driscoll and Heryudomo "Adaptive subsampling methods for RBF
% interpolation and collocation problems" 2007
% Diferente planteamiento a la hora de seleccionar los nuevos centros  
% a los existentes: se utiliza el valor del interpolante en ellos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Función a interpolar
f= @(x) 1./(1+25*x.^2);
%f= @(x) tanh(1000*(x-0.5));
%f= @(x) tanh(60*(x-0.01));

%% RBF que se usa
phi=@(r,epsilon) sqrt((epsilon*r).^2+1);

%% Intervalo de interpolación
a=-1;b=1;

%% Número inicial de centros
N=3;

%% Conjunto inicial de centros
x=linspace(a,b,N)';
 
%% Tolerancia para el refinamiento
thetar=1e-6;

%% Proceso de refinamiento
refine=true;
while any(refine)
    N=length(x);dx=diff(x);
    epsilon=0.75*min([Inf;1./dx],[1./dx;Inf]);
    y=x(1:N-1)+0.3*dx;y=[y;1];
    z=x(1:N-1)+0.5*dx;
    
    A=zeros(N);B=zeros(N);C=zeros(N-1,N);D=zeros(N-1,N);
    for j =1:N
        A(:,j)=phi(x-x(j),epsilon(j));
        B(:,j)=phi(y-y(j),epsilon(j));
        C(:,j)=phi(z-x(j),epsilon(j));
        D(:,j)=phi(z-y(j),epsilon(j));
    end
lambda1=A\f(x);lambda2=B\f(y);resid=abs(C*lambda1-D*lambda2);
refine=resid>thetar;x=sort([x;z(refine)]);
end

%% Resultados 
disp('Condicionamiento de la matriz')
disp(cond(A))
plot(x,0,'bo',x,A*lambda1,'k+',x,f(x),'ro')
disp('Número final de centros')
disp(length(x))

%% Malla de puntos para calcular el error
M=max(2*N,10^4);
 zz=linspace(a,b,M)';
 C=zeros(M,N);
 for j =1:N
       C(:,j)=phi(zz-x(j),epsilon(j));
 end
 disp('Error máximo sobre la malla')
 disp(max(abs(C*lambda1-f(zz))))

