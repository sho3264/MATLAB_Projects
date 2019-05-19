% Modal analysis
% Takes in Mass matrix M, Damping matrix C, Stiffness matrix K, External
% force F, initial conditions x0 and xdot0
% M,C, and K are square matricies of same size
% Try free vibration only
% Force vibration use ode45 to solve instead
function ff=modalode(M,C,K,x0,xdot0,F)
n=length(M);
% Make sure that all the input matrix have the same size
if ((n~=length(C))||(n~=length(K))||n~=length(F))
    ff=0;
else
    
    n=length(M);
    % Create symbolic variables w t
    % w is the frequencies
    % t is for time
    % Can symbolically plot it in 3D and use get frame better to do 3D plot
    syms w;
    
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
    % Change forces in x space to g space
    newF=U'*F;
    
    % Get transformed initial conditions
    % g0 is displacement at time t=0;
    g0=U'*M*x0;
    gdot0=U'*M*xdot0;
    
    % Create structured array to hold different time and displacements
    % First column represents time and second column represents
    % displacement corresponding to the time values in the same row in
    % the structured array
    % Each row in the structured array would be a different
    % mode/vibrational state of the system
    
    % Structured array in x space
    ff(1)=struct('time',[],'dispx',[]);
    
    
    % Create holder used to eventually combine all displacements
    gg(n)=struct('dispg',[]);
    
    % Have it run from 0 to 10 seconds
    tspan=[0 10];
    for pp=1:n
        % Create placeholders for initial conditions and all other
        % variables
        % Need to convert the values to doubles
        % If taken directly from matrix, the values are symbolic and will
        % make ode45 fail
        z0=[double(g0(pp,1)); double(gdot0(pp,1))];
        ww=double(wn(pp,1));
        psip=double(psi(pp,1));
        % Do not need to convert
        newFF=newF(pp,1);
        [tt,zz]=ode45(@(tt,zz) myfcn(tt,zz,ww,psip,newFF),tspan,z0);
        
        
        gg(pp).dispg=zz(:,1)';
        
        % Initialize values on first run only
        if pp==1
            % Cannot add empty arrays to other arrays
            hh=zeros(1,length(tt));
            ff(1).time=tt;
            holder=length(ff(1).time);
            vv=gg(pp).dispg;
        else
            % Make sure that the smallest matrix size is the size for all t
            % values
            % If current time matrix and previous time matrix are of the
            % same size do nothing
            if holder==length(tt)
                % No change needed
            else
                % Check if current or previous matrix size is smaller
                holder=min([holder,length(tt)]);
                ct=ff(1).time;
                ct1=length(ct);
                hh=zeros(n,holder);
                % If holder is smaller than last iteration of for loop, then
                % must change time part of structured matrix
                if holder==length(tt)
                    % Need to rewrite all previous matrix sizes to fit current
                    % one
                    % Temp variable to hold new created matricies
                    newhold=[];
                    for qq=1:pp-1
                        
                        % Another temp variable for old displacement
                        newhold1=ff(1).dispx(qq,:);
                        
                        % Modified displacement matrix to fit modified time
                        % variables
                        % Make sure first and last points are same
                        hh(:,1)=newhold1(:,1);
                        hh(:,holder)=newhold1(:,ct1);
                        
                        % Replace old time values with new time values
                        ff(1).time=tt;
                        vv=gg(pp).dispg;
                        % Must adjust previous held displacement matrix to have
                        % displacemetns corresponding to these new time values
                        for spec=2:holder-1
                            comp1=tt(spec,1);
                            for timecomp1=2:ct1
                                if comp1>=ct(timecomp1-1,1) && comp1<=ct(timecomp1,1)
                                    
                                    % Check which time value the comparison
                                    % number is closer to
                                    c1=comp1-ct(timecomp1-1,1);
                                    c2=ct(timecomp1,1)-comp1;
                                    
                                    if c1==c2
                                        % Comparison number is half way between
                                        % 2 time points in old time matrix
                                        hh(:,spec)=(newhold1(:,timecomp1)+newhold1(:,timecomp1-1))/2;
                                    elseif c1<c2
                                        % Comparison number is closer to
                                        % smaller time value
                                        hh(:,spec)=newhold1(:,timecomp1-1);
                                    else
                                        % Comparison number is closer to larger
                                        % time value
                                        hh(:,spec)=newhold1(:,timecomp1);
                                    end
                                    % Replace old displacement values with new
                                    % ones
                                else
                                end
                            end
                        end
                        hh=hh(:,1:holder);
                        newhold=[newhold;hh];
                    end
                    
                    ff(1).dispx=newhold(1:(pp-1),:);
                else
                    vv=zeros(1,holder);
                    % Make sure first and last points are same
                    vv(:,1)=gg(:,1);
                    vv(:,holder)=gg(:,length(gg));
                    
                    % Must change new displacement matrix to fit old time values
                    for spec2=2:holder-1
                        comp2=ct(spec2,1);
                        
                        for timecomp2=2:length(zz)
                            % Try to have new time values come close as
                            % possible to old ones
                            if comp2>=tt(timecomp2-1,1) && comp2<=tt(timecomp2,1)
                                
                                % Check which time value the comparison
                                % number is closer to
                                c3=comp2-tt(timecomp2-1,1);
                                c4=tt(timecomp2,1)-comp2;
                                if c3==c4
                                    % Comparison number is half way between
                                    % 2 time points in old time matrix
                                    % Take midpoint of displacement
                                    vv(:,spec2)=(newhold1(:,timecomp2)+newhold1(:,timecomp2-1))/2;
                                elseif c3<c4
                                    % Comparison number is closer to
                                    % smaller time value
                                    vv(:,spec2)=newhold1(:,timecomp2-1);
                                else
                                    % Comparison number is closer to larger
                                    % time value
                                    vv(:,spec2)=newhold1(:,timecomp2);
                                end
                            else
                            end
                            
                        end
                        
                    end
                end
                
            end
        end
        
        
        % Sum up all sets of points from the contribution of all
        % modes/states
        ff(1).dispx=[ff(1).dispx;vv]
        
    end
    % Change back to x space from g space
    ff(1).dispx=(U*ff(1).dispx)'
    
end
end



