%Ordering points
%Takes in array and reorders rows based on first element of each row
function l=order(a)
%Find length of array
g=length(a);
    for i=1:g-1
        
        for k=i:g-1
            %Use placeholder variable to swap with other
            if a(k,1)>a(k+1,1)
                %placeholder variable m
                m=a(k,:);
                a(k,:)=a(k+1,:);
                a(k+1,:)=m;
            end
           
        end

    end
    
l=a;
    
end
