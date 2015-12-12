


%%%%% the function solves the problem min(f'x) sub to ax<=b.
% x returns optimal point in tabulaeu representation. 
% fval returns the min(f'x)
% it returns the iteration.
% op=1 means optimal solution exists. otherwise doesn't
function [x,fval,it,op]=singl(f,a,b) 
[m,n]=size(a);  
c=[a eye(m) b; f' zeros(1,m+1)]; 
fval=0;  
x=zeros(m+n,1);

op=1;
it=0; 
e=zeros(1,m);
lie=find(f<0); 
l=length(lie); 
while(l>0)   
    
    for j=1:l    
        d=find(c(:,lie(j))); 
        
        d_l=length(d);  
        if d_l>0         
            for i=1:m    
                if c(i,lie(j))>0  
                    e(i)=c(i,end)/c(i,lie(j));  
                else
                    e(i)=inf;       
                end
            end
            [g,h]=min(e); 
            for w=1:m+1  
                if w==h    
                    c(w,:)=c(w,:)/c(h,lie(j)); 
                else
                    c(w,:)=c(w,:)-c(h,:)*c(w,lie(j))/c(h,lie(j));   
                end
            end
            c
            it=it+1;   
        else
            op=0;    
        end
    end
    lie=find(c(end,:)<0);     
    l=length(lie);
    
    
    
    
    
end
for i=1:(m+n)   
    ix=find(c(:,i));   
    if(length(ix)==1)&(ix<=m)&(c(ix,i)==1)     
        x(i)=c(ix,end);
    else
        x(i)=0   ;
    end
end
fval=-c(end,end);
end