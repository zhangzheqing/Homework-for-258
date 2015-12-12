A=[1,0,0;0,1,0;1,1,1;0,1,3,;-1,0,0;0,-1,0;0,0,-1;1,-1,1;1,-1,-1;-1,1,-1.1];
b=[200;300;400;600;0;0;0;250;150;250];
f=-[4,-5,1]';
c=[4,-5,1];
f=[-1,-1,1]; % (0,0,200)
u=cross(c,f);
u=u/norm(u);
figure
polytope2D= shadow(A,b,c,f)
innerp=polytope2D(2:end,1:end)*c';
max(innerp)
min(innerp)
innerp
plotpoly3D(A,b)
count=0;

axis([-100 400 -100 400 -100 400])
for i= (-300):0.1:100
    count=count+1;
    Q(count,:)= [-62.5931  221.7586  -95.3800]+i*u;
    P(count,:)= [159.0164  -51.2295  -42.2131]+i*u;
    
    
end
plot3(Q(:,1),Q(:,2),Q(:,3),'g')
plot3(P(:,1),P(:,2),P(:,3),'b')

