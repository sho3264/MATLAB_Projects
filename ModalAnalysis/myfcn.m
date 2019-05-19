% Turn to first order diff eq
% p,x are needed to use ode45
% w is the frequency
% psi is the damping ratio
% f is the applied external force
function dz=myfcn(p,x,w,psi,f)
% Check if f is a symbolic function or a double
% If it is a symbolic function, must change it to a double
if strcmp(class(f),'sym')==1
    syms t;
    % Assuming input function f is a function of time t
    tspace=linspace(0,10,25);
    % Changing f into a double
    gg=double(subs(f,t,tspace));
    % 1-D interpolates the function f at points tspace
    ff=interp1(tspace,gg,p);
else
    ff=f;
end
dz=zeros(2,1);
dz(1)=x(2);
dz(2)=-w^2*x(1)-2*w*psi*x(2)+ff;
