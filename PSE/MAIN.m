clear all
%%% This is a  numerical integration program based on third-order Runge-Kutta methods and central difference.
%%% It works well with functions without non-analytic points. 
%%% However, if there are non-analytic points, the simulation will be somewhat ineffective.
%%% I am currently looking for a good method.

% *define the global variables and others we need

global U V alpha beta numPoints_phi_1 numPoints_theta_1 theta_2 dphi_1 phi_1 phi_2
theta_2 = pi/4; 
phi_2 = 0;

% the number of theta_1
numPoints_theta_1 = 200; 

% the number of phi_1
numPoints_phi_1 = 200; 

% the total time is T=  dt*T_num
dt = 0.001;
T_num = 600; 

dtheta_1 = pi/numPoints_theta_1; 
dphi_1 = 2*pi/numPoints_phi_1;
alpha = dt/(dphi_1);
beta = dt/(dtheta_1);
theta_1 = linspace(0,pi,numPoints_theta_1+1);
phi_1 = linspace(0,2*pi,numPoints_phi_1+1);

%Set up the grid
[Phi_1,Theta_1] = meshgrid(phi_1,theta_1); 
R = 1;

% change to the Cartesian coordinate system
X = R*sin(Theta_1).*cos(Phi_1);
Y = R*sin(Theta_1).*sin(Phi_1);
Z = R*cos(Theta_1);

%Compute the vectors
U = sin(theta_2)*cot(Theta_1).*cos(Phi_1-phi_2)-cos(theta_2);
V = sin(theta_2).*sin(Phi_1-phi_2);

% define the phase density at T=0 
pd0 = zeros([numPoints_theta_1+1,numPoints_phi_1+1]);
for x = [1:numPoints_theta_1+1]
    for y  = [1:numPoints_phi_1+1]
        pd0(x,y) = sin(theta_1(x)/2+pi/2).^2*cos(phi_1(y)/2+pi/2)^2; % This initial condition has a non-analytic points, So It doesn't perform well
    end
end
pd0 = X.^2.*Y.^2.*Z.^2; % change to another form

% Three order Runge-Kutta methods
pd = pd0;
for Num = [1:T_num]
    pd1 = Evolution(pd);
    pd2 = 3/4*pd+1/4*Evolution(pd1);
    pd = 1/3*pd+2/3*Evolution(pd2);
end

% Plot the result 
figure
set(gcf,'unit','centimeters','position',[3 3 25 15]);
t=tiledlayout(1,2,'TileSpacing','Compact');  

nexttile
FS = 20;
surf(X,Y,Z,pd0);
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
yticks([-1,1])
yticklabels({'-1','1'})
ylabel('$$y$$','Interpreter','latex')
xticks([-1,1])
xlabel('$$x$$','Interpreter','latex')
xticklabels({'-1','1'})
zticks([-1,1])
zticklabels({'-1','1'})
zlabel('$$z$$','Interpreter','latex')
shading interp
set(gca,'TickLabelInterpreter','latex'); 
title(['$$time = 0$$'],'Interpreter','latex','FontSize',FS+2)
set(gca,'FontSize',FS);
set(gca,'YAxisLocation','right'); 

nexttile
surf(X,Y,Z,pd);
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
yticks([-1,1])
yticklabels({'-1','1'})
ylabel('$$y$$','Interpreter','latex')
xticks([-1,1])
xlabel('$$x$$','Interpreter','latex')
xticklabels({'-1','1'})
zticks([-1,1])
zticklabels({'-1','1'})
zlabel('$$z$$','Interpreter','latex')
shading interp
set(gca,'TickLabelInterpreter','latex'); 
title(['$$time =0.6$$'],'Interpreter','latex','FontSize',FS+2)
set(gca,'FontSize',FS);
set(gca,'YAxisLocation','right'); 

% C = colorbar('position',[0.1,0.1,0.8,0.01],'orientation','horizontal')
% set(C,'Ytick',0:0.02:0.04)
% 
% set(C,'TickLabelInterpreter','latex');
% caxis([0,0.04])
% 
% set(gca,'FontSize',FS);
% set(gca,'YAxisLocation','right'); 
f = gcf;
exportgraphics(f,'fig1.pdf')