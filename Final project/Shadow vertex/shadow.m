
function polytope2D= shadow(A,b,c,f)
% projection from Ax <= b to span(c,f) 
% polytope2D returns the vertex of the projection. The projection of a 
% polytope on 2D plane is a polygon.
u=cross(c,f);
vertex=Vertex(A,b);
l=length(vertex(:,1));

polytope2D=[pi,pi,pi];
for i = 1:(l-1)
    for j=(i+1):l
        ver1=vertex(i,:);
        ver2=vertex(j,:);
        vmid=1/2 * (ver1+ver2);
        error=b-A*vmid';
        inter_plane=[];
        count=0;
        
        for kk=1:length(error)
            if error(kk)<10^(-11);
                count=count+1;
                inter_plane(count,:)=A(kk,:);
            end
        end
            
        if count ==2
            cone= (u*inter_plane(1,:)')*(u*inter_plane(2,:)');
            if cone <= 0   % ver1 and ver2 appear in the shadow.
                count1=0;
                count2=0;
                %%%%%% calculate the projection of ver1 ver2
                matr=[c*c',c*f';c*f',f*f'];
                br=[c*ver1';f*ver1'];
                lambda=matr\br;
                ver1plane=c*lambda(1)+f*lambda(2);

                br=[c*ver2';f*ver2'];
                lambda=matr\br;
                ver2plane=c*lambda(1)+f*lambda(2);
                
               
                for kkk=1:length(polytope2D(:,1))
                     
                      diff1=norm(polytope2D(kkk,:)-ver1plane);
                      if diff1 <10^(-10)
                          count1=count1+1;
              
                      end
                      
                      
                end
                
                if count1 ==0
                    polytope2D =[ polytope2D; ver1plane];
                end
                    
                for kkk=1:length(polytope2D(:,1))
                    
                      diff2=norm(polytope2D(kkk,:)-ver2plane);
                      if diff2 <10^(-10)
                          count2=count2+1;
                          
                      end
                end
                
                if count2 ==0
                    polytope2D =[ polytope2D; ver2plane];
                end
                
                
                
                
                
            end
            
        end
    end
    
end



plot3(polytope2D(2:end,1),polytope2D(2:end,2),polytope2D(2:end,3),'o')
xlabel('x1')
ylabel('x2')
zlabel('x3')



end