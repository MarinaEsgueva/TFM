function front=u2_cond_frontera(x,y)
%Condiciones de frontera para la  EDP 2
%Cuando x=1  0 x=0 devuelve u(x,y)
%Cuando y=0 o y=1 devuelve u_y(x,y)
if x==0
    front=1;
elseif x==1
    front=0.1;
elseif y==0|y==1
    front=0;
end