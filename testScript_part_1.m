clear all
close all
warning off
format short g
%% *****PART 1*****
gamma = 0.1:0.001:0.2;
Y=zeros(4,length(gamma));
options = optimoptions('fsolve','Display','none');

P_plot=zeros(length(gamma),8);
Z_plot=zeros(length(gamma),8);
N_plot=zeros(length(gamma),8);
D_plot=zeros(length(gamma),8);

X(:,1) = [0.85; 0.05; 0.05; 0.05];
X(:,2) = [0.05; 0.85; 0.05; 0.05];
X(:,3) = [0.05; 0.05; 0.85; 0.05];
X(:,4) = [0.05; 0.05; 0.05; 0.85];
X(:,5) = [0.55; 0.15; 0.15; 0.15];
X(:,6) = [0.15; 0.55; 0.15; 0.15];
X(:,7) = [0.15; 0.15; 0.55; 0.15];
X(:,8) = [0.15; 0.15; 0.15; 0.55];

for i = 1:length(X(1,:))
    A = X(:,i);
    for j = 1:length(gamma)
        setGlobalx(gamma(j));
        Y(:,j) = fsolve(@func, A, options);
    end
    P_plot(:,i) = (Y(1,:)).';
    Z_plot(:,i) = (Y(2,:)).';
    N_plot(:,i) = (Y(3,:)).';
    D_plot(:,i) = (Y(4,:)).';
end
%%
% finding \gamma^*
disp('*********************************************************')
disp('Compare to bifurcations to verify these values of gamma^*')
P_plot(38,:).'
Z_plot(38,:).'
N_plot(38,:).'
D_plot(38,:).'
%%
% bifurcations
figure
hold on
title('phytoplankton')
ylabel('{grams of carbon}/{meter^3}');
xlabel('\Delta \gamma')
for i=1:length(X(1,:))
    plot(gamma,P_plot(:,i),'-b');
end
hold off
figure
hold on
title('zooplankton')
ylabel('{grams of carbon}/{meter^3}');
xlabel('\Delta \gamma')
for i=1:length(X(1,:))
    plot(gamma,Z_plot(:,i),'-b');
end
hold off
figure
hold on
title('nutrients')
ylabel('{grams of carbon}/{meter^3}');
xlabel('\Delta \gamma')
for i=1:length(X(1,:))
    plot(gamma,N_plot(:,i),'-b');
end
hold off
figure
hold on
title('detritus')
ylabel('{grams of carbon}/{meter^3}');
xlabel('\Delta \gamma')
for i=1:length(X(1,:))
    plot(gamma,D_plot(:,i),'-b');
end
hold off
%%
% **periodic behavior**
% finding min--max values on the limit cycle
x0=[0.15; 0.06; 0.65; 0.14];
gamma_s=0.12:0.0001:0.2;
time = 0:1000;
min_max_P=zeros(length(gamma_s), 2);
min_max_Z=zeros(length(gamma_s), 2);
min_max_N=zeros(length(gamma_s), 2);
min_max_D=zeros(length(gamma_s), 2);

%figure;
%hold on
for i=1:length(gamma_s)
    [t, y] = ode45(@RHS_eqs, time, [x0; gamma_s(i)]);
    %{
    subplot(5,1,i)
    plot(t,y,'-');
    title(strcat('for \gamma =',num2str(gamma_s(i))));
    xlabel('time');
    ylabel('{gC}/{m^3}');
    L=legend('phytoplankton','zooplankton','nutrients','detritus','\gamma');
    set(L, 'FontSize', 6)
    %}
    min_max_P(i,:)= [min(y(:,1)), max(y(:,1))];
    min_max_Z(i,:)= [min(y(:,2)), max(y(:,2))];
    min_max_N(i,:)= [min(y(:,3)), max(y(:,3))];
    min_max_D(i,:)= [min(y(:,4)), max(y(:,4))];
end
hold off
%%
disp('*********************************************************')
disp('*******pulling min/max values from ode45 function********')

figure
axis([0.12 0.201 -0.1 1.01])
hold on
plot(gamma_s, min_max_P(:,1),'-r', 'LineWidth', 1.5);
plot(gamma_s, min_max_Z(:,1),'-b', 'LineWidth', 1.5);
plot(gamma_s, min_max_N(:,1),'-g', 'LineWidth', 1.5);
plot(gamma_s, min_max_D(:,1),'-c', 'LineWidth', 1.5);
plot(gamma_s, min_max_P(:,2),'-r', 'LineWidth', 1.5);
plot(gamma_s, min_max_Z(:,2),'-b', 'LineWidth', 1.5);
plot(gamma_s, min_max_N(:,2),'-g', 'LineWidth', 1.5);
plot(gamma_s, min_max_D(:,2),'-c', 'LineWidth', 1.5);
title('min/max values')
legend('phytoplankton','zooplankton','nutrients','detritus');
xlabel('\gamma')
ylabel('gC/m^3')
hold off

%}
%%
% ***Some behavior testing***
%{
figure;
hold on
[t, y] = ode45(@RHS_eqs, 0:10000, [x0; 0.134]);
plot(t,y,'-b');
hold off
%}
%{
gamma_st = 0.1365:0.0001:0.1372;
for i=1:10
    hold on
    [t, y] = ode45(@RHS_eqs, 0:5000, [x0; gamma_st]);
    plot(t,y,'-r');
    hold off
end
%}
%% ***Functions***
function F = func(y)
k_N=0.2;
k_P=0.1;
lambda_P=0.1;
I=0.6;
eps=0.3;
beta=0.4;
lambda_Z=0.1;
nu=0.5;
phi=0.1;
gamma = getGlobalx;

P = y(1);
Z = y(2);
N = y(3);
D = y(4);

dP=(N/(k_N+N))*gamma*P-lambda_P*P-(P/(k_P+P))*I*Z;
dZ=eps*(P/(k_P+P))*I*Z-lambda_Z*Z;
dN=phi*D+beta*(P/(k_P+P))*I.*Z-(N/(k_N+N))*gamma*P+nu*lambda_Z*Z;
dD=-phi*D+lambda_P*P+(1-eps-beta)*(P/(k_P+P))*I*Z+(1-nu)*lambda_Z*Z;
d_sum=P+Z+N+D-1;

F = [dP; dZ; dN; dD; d_sum];
end

function setGlobalx(val)
global x
x = val;
end

function r = getGlobalx
global x
r = x;
end

