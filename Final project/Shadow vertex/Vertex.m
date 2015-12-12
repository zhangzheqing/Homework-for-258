function vertex = Vertex(A,b)

% "vertex" returns the matrix recording the coordinates of the vertex. 
l=length(A(:,1));
vertex=[];
num=0;

for i = 1:(l-2) 
    for j=(i+1):(l-1)
        for k =(j+1):l
            B=[A(i,:);A(j,:);A(k,:)];
            bk=[b(i),b(j),b(k)]';
            vert=B\bk;
            diff=A*vert-b;
            m=max(diff);
            if m <= 0 
                num=num+1;
                vertex(num,:)=vert';
                
            end
        end
    end
end


end


         
            
            
        









