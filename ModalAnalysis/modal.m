% Modal analysis
% Takes in Mass matrix M, Damping matrix C, Stiffness matrix K, External
% force F, initial conditions x0 and xdot0
% M,C, and K are square matricies of same size
% Try free vibration only
% Force vibration use ode45 to solve instead
function ff=modal(M,C,K,x0,xdot0)
n=length(M);
% Make sure that all the input matrix have the same size
if (n~=length(C))||(n~=length(K))
    ff=0;
else
    
    
    % Create symbolic variables w t
    % w is the frequencies
    % t is for time
    % Can symbolically plot it in 3D and use get frame better to do 3D plot
    syms w t;
    
    %Solve for frequencies
    ll=-w^2*M+K;
    % Reorder the frequencies from smallest to largest
    f=order(solve(det(ll)==0,w));
    % Natrual frequencies of the system
    wn=[];
    X=[];
    
    for i=1:length(f)
        % Make sure to get only positive frequencies
        if f(i,1)>0
            wn=[wn;f(i,1)];
        else
        end
    end
    % Make sure to organize all frequencies and consider possibility of
    % frequencies equal to 0
    wn=order(wn);
    % See if only single mass system
    % eig of single mass system yields 0
    if n~=1
        for ii=1:n
            % Create temp variable for eig functions and eig variables
            % Create matrix Ef which is eigenfunctions corresponding to
            % eigenvalues in matrix Ev
            [Ef,~]=eig(subs(ll,w,wn(ii,1)));
            %Create eigenvector matrix
            X=[X Ef(:,1)];
        end
        
        % Create alpha matrix to multiply with eigenvector matrix to form mass
        % Normalized eigenvectors
        alpha=[];
        for jj=1:n
            % Temp variable to hold quantity and append to
            a11=1/(X(:,jj)'*M*X(:,jj))^.5;
            alpha=[alpha;a11];
        end
        
        % Mass normalized eigenvectors
        U=[];
        for count=1:n
            
            U=[U alpha(count,1)*X(:,count)];
            
        end
        
    else
        
        U=1;
    end
    % Uncouple equations
    % Create psi matrix if there is damping
    if C==zeros(n)
        psi=zeros(n,1);
    else
        psi=[];
        
        for count2=1:n
            % Another Temp variable used for placeholder
            a12=U'*C*U;
            psi=[psi;a12(count2,count2)/(2*wn(count2,1))];
        end
        
    end
    % Create damped frequencies
    wd=zeros(n,1);
    
    for count3=1:n
        
        wd(count3,1)=wn(count3,1)*(1-psi(count3,1)^2)^.5;
        
    end
    
    % Get transformed initial conditions
    % g0 is displacement at time t=0;
    g0=U'*M*x0;
    gdot0=U'*M*xdot0;
    
    % Assume sine + cosine solution and solve for coefficients
    % A is associated with the cosine term and B with the sine term
    A=zeros(n,1);
    B=zeros(n,1);
    
    for count4=1:n
        
        A(count4,1)=g0(count4,1);
        B(count4,1)=(gdot0(count4,1)+psi(count4,1)*wn(count4,1)*g0(count4,1))/wd(count4,1);
        
    end
    % Combine everything and print results as column vector
    ff=[];
    for pp=1:n
        % Another Temp variable placeholder
        % Basics form of solution to Free vibration with damping
        % If no damping, psi=0, therfore exp term goes to 1
        gg=exp(-psi(pp,1)*wn(pp,1)*t)*(A(pp,1)*cos(wd(pp,1)*t)+B(pp,1)*sin(wd(pp,1)*t));
        
        ff=[ff; gg];
    end
    ff=U*ff;
end
end


