% demo for LSLDG clustering

clear all;

rng('shuffle');

op.dim=2; % data dimension 
op.samples=1000; % number of samples
op.bnum=min(100,op.samples); % number of basis functions
op.cvfold=5; % N-fold cross validation
op.sigma_list=logspace(-1,1,6); % list of sigma used in cross validation
op.lambda_list=logspace(-2,1,6); % list of lambda used in cross validation
op.maxiter=100; % maximum number of iteration in clustering update
op.tol=1e-10; % Stopping criteria in clustering

% Data for a mixture of three Gaussians
c=2; s=1;
myu = [0 c zeros(1,op.dim-2);-c -c zeros(1,op.dim-2);c -c zeros(1,op.dim-2)]; 
Smat = cat(3,eye(op.dim),eye(op.dim),eye(op.dim)); 
p = [0.4,0.3,0.3]; GMMobj = gmdistribution(myu,Smat,p);
X=random(GMMobj,op.samples)';

figure(1); subplot(121); scatter(X(1,:),X(2,:),'ko'); 
axis([-5 5 -5 5]); axis square; 
title('Input Data','FontSize',15); hold on;
plot(0,c,'r*'); plot(-c,-c,'b*'); plot(c,-c,'g*'); 
set(gca,'FontSize',15); hold off;

% LSLDG Clustering
op.disp=1; % Display the transition of data points during clustering updates
Y=LSLDGClustMatrix(X,op); % Matrix-form update


figure(1); subplot(122); scatter(Y(1,:),Y(2,:),'ko'); 
axis([-5 5 -5 5]); axis square; 
title('LSLDG Clustering','FontSize',15); hold on;
plot(0,c,'r*'); plot(-c,-c,'b*'); plot(c,-c,'g*'); 
set(gca,'FontSize',15); hold off;