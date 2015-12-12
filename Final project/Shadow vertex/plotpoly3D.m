
%A=[1,0,0;0,1,0;1,1,1;0,1,3,;-1,0,0;0,-1,0;0,0,-1;1,-1,1;1,-1,-1;-1,1,-1.1];
%b=[200;300;400;600;0;0;0;250;150;250];
function plotpoly3D=plotpoly3D(A,b)
vertex=Vertex(A,b);
l=length(vertex(:,1));
num_line=0;
for i = 1:(l-1)
    for j=(i+1):l
        ver1=vertex(i,:);
        ver2=vertex(j,:);
        vmid=1/2 * (ver1+ver2);
        error=b-A*vmid';
        count=0;
        
        for kk=1:length(error)
            if error(kk)<10^(-17);
                count=count+1;
            end
        end
            
        if count >=2
            num_line=num_line+1;
            n=0;
            line12=[];
            for lambda=0:0.1:1
                n=n+1;
                line12(n,:)= ver1*lambda+ver2*(1-lambda);
            end
            plot3(line12(:,1),line12(:,2),line12(:,3),'r');
            hold on 
        end
    end
    
end

xlabel('x1')
ylabel('x2')
zlabel('x3')
end