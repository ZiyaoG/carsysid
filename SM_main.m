% This script aims at taking a csv file and do set membership estimation
% based on the data. The data is read from a experimental car. 
clear all;
close all;
yalmip('clear')
%% read data and parameters
data1=readmatrix("log_state.csv");
vx=data1(:,1);
vy=data1(:,2);
psi_dt=data1(:,3);
theta=data1(:,4);
phi=data1(:,5);
delta=data1(:,6);

% some known parameters
g=9.8;
f=30;
dt=1/f;

% prior knowledge
A0=[-4 -12; 1 -5];
B0=[24; 15];
w_bound=1e10*[1;1]; % bound on the "lumped-up" disturbances/noises

% parameters
N=500; % buffer size
N_start=10;
%% LSE update A_LSE, B_LSE




%% set membership update A_SM, B_SM
d=[-g*sin(phi).*cos(theta),zeros(size(-g*sin(phi).*cos(theta)))];
d=d(N_start+1:N_start+N,:);
z=[vx psi_dt];
zz=z(N_start+2:N_start+N+1,:);
z=z(N_start+1:N_start+N,:);
%% debug
% figure(1)
% subplot(2,1,1)
% plot(1:N,zz(:,1)-z(:,1));
% subplot(2,1,2)
% plot(1:N,zz(:,2)-z(:,2));
% 
% return;
%%
u=delta(N_start+1:N_start+N);
vx=vx(N_start+1:N_start+N,:);

% unfalsified sets
A_ij=sdpvar(2,2); % NOTE: THIS IS THE A_ij AFTER EXCLUDING V_x
B_est=sdpvar(2,1);
constraints=[-A_ij(1,1)<=-1e-7, -A_ij(2,2)<=-1e-7, -B_est(1,1)<=-1e-7, -B_est(2,1)<=-1e-7];
for i=1:N
    vx_temp=vx(i,:);
    A_est=[-A_ij(1,1)/vx_temp -vx_temp-A_ij(1,2)/vx_temp; -A_ij(2,1)/vx_temp  -A_ij(2,2)/vx_temp];
    cons=zz(i,:)'-z(i,:)'-(d(i,:)'+A_est*z(i,:)'+B_est*u(i,1))*dt<=w_bound;
    constraints=[constraints,cons];
end

A_distance=A0-A_est;
B_distance=B0-B_est;
% Objective = max(abs(Distance(:))); % minimize the infty norm of difference
Objective = norm(A_distance)+norm(B_distance); % minimize 2-norm
options = sdpsettings('verbose', 0, 'solver', 'mosek', ...
                              'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 100, ... % 300 seconds time limit
                              'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-4, ... % Adjust relative gap tolerance
                              'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-4); % Adjust co-tolerance
        sol = optimize(constraints, Objective, options);
        % sol = optimize(Constraints, Objective);


        if sol.problem == 0
            A_est_value = double([-A_ij(1,1)/vx_temp -vx_temp-A_ij(1,2)/vx_temp; -A_ij(2,1)/vx_temp  -A_ij(2,2)/vx_temp]);
            B_est_value = double(B_est);
            disp('-------------------- Solution Found------------------------')
            % Print A_est
            fprintf('A_est =\n');
            for i = 1:size(A_est_value, 1)
                fprintf('   ');
                fprintf('%.1f  ', A_est_value(i, :)); % Print each row
                fprintf('\n');
            end
            
            % Print B_est
            fprintf('B_est =\n');
            for i = 1:size(B_est_value, 1)
                fprintf('   ');
                fprintf('%.1f  ', B_est_value(i, :)); % Print each row
                fprintf('\n');
            end
        else
            disp('No solution found');
            return;
        end
%% plot
% B_polytope=Polyhedron()
% plot(-A_ij(1,1)/vx_temp -vx_temp-A_ij(1,2)/vx_temp)