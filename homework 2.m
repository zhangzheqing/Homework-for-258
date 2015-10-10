%%% homework 2 problem 3  by Zheqing Zhang

%%% I am using gradient descent method with backtracking linesearch to
%%% solve min(- (a u + b)+log(1+exp(a u + b))). Since the logistic function
%%% is convex, we can use first order gradient descent method to solve.



M=csvread('binary.csv',1,0);

m=length(M(:,2));

%normalization
for i=2:3
    M(:,i)=(M(:,i)-min(M(:,i)))/(max(M(:,i))-min(M(:,i)));
end
% start point u_intial=[2,2] beta_intial=2.
a1=2;
a2=2;
beta=2;

%backtracking parameters.
alpha=0.2;
beta0=0.5;



%count the steps taken. 
in=1; 


%save the data of u and beta at each step.
X=[a1;a2;beta];




% calculate the grad at (a1 a2 beta) 
grad=zeros(3,1);

for i = 1 : m
    
    grad(1) = grad(1) - M(i,1)*M(i,2)+(1-M(i,1))*...
        exp(a1*M(i,2)+a2*M(i,3)+beta)*M(i,2)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
    
    grad(2) = grad(2) - M(i,1)*M(i,3)+(1-M(i,1))*...
        exp(a1*M(i,2)+a2*M(i,3)+beta)*M(i,3)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
    
    grad(3) = grad(3) - M(i,1)+(1-M(i,1))* ...
        exp(a1*M(i,2)+a2*M(i,3)+beta)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
end

dx= -grad/norm(grad,2); % the steepest descent direction.

norm_grad=norm(grad,2);
G=[norm_grad];

while norm_grad >  0.00001 %%%%%%%% the stop criteria 
    
  

in=in+1;    % number of iteration increased by 1.



t=1;

    L_left=1;
    L_right=0;
    count=0;
    while   L_left > L_right    % f(x+t dx) > f(x)+ alpha* t* grad(f)* dx  
    count=count+1;
   
        % calculate f(x+t dx)
    L_left=0;
    L_right=0;
    vector=X(:,in-1)+dx*t;

    for i= 1: m
        L_left= L_left-M(i,1)*(vector(1)*M(i,2)+vector(2)*M(i,3)+vector(3))+(1-M(i,1))*...
        log(1+exp(vector(1)*M(i,2)+vector(2)*M(i,3)+vector(3)));
    
        L_right= L_right-M(i,1)*(X(1,in-1)*M(i,2)+X(2,in-1)*M(i,3)+X(3,in-1))+(1-M(i,1))*...
        log(1+exp(X(1,in-1)*M(i,2)+X(2,in-1)*M(i,3)+X(3,in-1)));
    end
    L_right=L_right-alpha*t*norm(grad,2);
    t=t*beta0;

    end
X(:,in)=X(:,in-1)+dx*t;

a1=X(1,in);
a2=X(2,in);
beta=X(3,in);

grad=zeros(3,1);
for i = 1 : m
    
    grad(1) = grad(1) - M(i,1)*M(i,2)+(1-M(i,1))*...
        exp(a1*M(i,2)+a2*M(i,3)+beta)*M(i,2)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
    
    grad(2) = grad(2) - M(i,1)*M(i,3)+(1-M(i,1))*...
        exp(a1*M(i,2)+a2*M(i,3)+beta)*M(i,3)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
    
    grad(3) = grad(3) - M(i,1)+(1-M(i,1))* ...
        exp(a1*M(i,2)+a2*M(i,3)+beta)/(1+exp(a1*M(i,2)+a2*M(i,3)+beta));
end


dx= -grad/norm(grad,2);
norm_grad=norm(grad,2);


G=[G,norm_grad];
T=[T,t];
end



a1=X(1,in);
a2=X(2,in);
beta=X(3,in);
x=0:0.01:1;
y=0:0.01:1;

p = @(u1,u2) exp(a1*u1+a2*u2+beta)/(1+exp(a1*u1+a2*u2+beta)); 

P=zeros(length(x),length(y));

for i=1:length(x)
    for j=1:length(y)
        P(i,j)=p(x(i),y(j));
    end
end
fprintf('a1=%f\na2=%f\nbeta=%f\n',a1,a2,beta)

figure(1)
mesh(y,x,P);
xlabel('GRE')
ylabel('GPA')
zlabel('admission')
title('logistic regression ')

figure(2)
M=csvread('binary.csv',1,0);
for i=1:m
    if M(i,1)>0.5
        plot(M(i,2),M(i,3),'g+')
        hold on
    else
        plot(M(i,2),M(i,3),'ro')
        hold on
    end
end

  xlabel('gre')
  ylabel('gpa')
  title('scatter plot of raw data')
  
  figure(3)
  
  
  x=0:0.01:1;
  y=[];
  for i=1:length(x)
      y(i)=(1/2-beta-a1*x(i))/a2
  end
  plot(x,y,':')
  
  xlabel('normalized gre')
  ylabel('normalized gpa')
  title('decision boundary')
    
    
    
    

















