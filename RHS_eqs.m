function output=RHS_eqs(t,X) 
%Note, t doesn't show up explicitly anywhere in our problem, but ODE45 needs to input both
%indepdend and dependent variables for considering the most general case.


%Book keeping for efficiency
M=length(X)/5;
P=X(1:M); 
Z=X(M+1:2*M); 
N=X(2*M+1:3*M); 
D=X(3*M+1:4*M); 
gamma=X(4*M+1:5*M);

%Other non-gamma parameters in the model
k_N=0.2; 
k_P=0.1;
lambda_P=0.1; 
I=0.6;
eps=0.3; 
beta=0.4; 
lambda_Z=0.1; 
nu=0.5; 
phi=0.1;

%Modulation function
f = @(x,k) x./(k+x);

%This is the meat of this function
P_RHS=f(N,k_N).*gamma.*P-lambda_P*P-f(P,k_P).*I.*Z;
Z_RHS=eps*f(P,k_P).*I.*Z-lambda_Z.*Z;
N_RHS=phi*D+beta*f(P,k_P).*I.*Z-f(N,k_N).*gamma.*P+nu*lambda_Z.*Z;
D_RHS=-phi*D+lambda_P*P+(1-eps-beta)*f(P,k_P).*I.*Z+(1-nu)*lambda_Z*Z;
%%
% More book keeping for output formating
%%
temp(1:M)=P_RHS; 
temp(M+1:2*M)=Z_RHS; 
temp(2*M+1:3*M)=N_RHS; 
temp(3*M+1:4*M)=D_RHS; 
temp(4*M+1:5*M)=zeros(1,M); %Note, d gamma /dt= 0
temp=temp';

%Right hand side evaluations above are our output
output=temp;


