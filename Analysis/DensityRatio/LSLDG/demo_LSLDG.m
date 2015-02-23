clear all;

rng('shuffle');

op.dim=1; % data dimesnion 
op.samples=1000; % number of samples
op.bnum=min(100,op.samples); % number of basis functions
op.cvfold=5; % N-fold cross validation
op.sigma_list=logspace(-1,1,6); % list of sigma used in cross validation
op.lambda_list=logspace(-2,1,6); % list of lambda used in cross validation

% Gaussian data
X=randn(op.dim,op.samples);

% log-density gradient estimation
g=LSLDG(X,op);

[xx,ind]=sort(X,'descend');
figure(1); subplot(121);
plot(xx,g(ind),'*b'); 
axis([-3 3 -3 3]); hold on;

xxt=linspace(-5,5,100000);
figure(1); subplot(121);
plot(xxt,-xxt,'-r','LineWidth',1.5);  
title('Gaussian Data','FontSize',15);
axis([-3 3 -3 3]); legend('LSLDG','True'); 
axis square; set(gca,'FontSize',15); hold off;



% data from a mixture of two Gaussians
c=2; MU=[c;-c]; SIGMA=cat(3,1,1); p1=0.5; p2=1-p1; p = [p1,p2]; 
obj = gmdistribution(MU,SIGMA,p); X = random(obj,op.samples)';

% log-density gradient estimation
g=LSLDG(X,op);

[xx,ind]=sort(X,'descend');
figure(1); subplot(122);
plot(xx,g(ind),'*b'); 
axis([-5 5 -3 3]); hold on;

% true log-density gradient
xxt=linspace(-5,5,100000);
pdfFig=p1*exp(-0.5*sum((xxt-c).^2,1))/sqrt(2*pi)^op.dim...
    +p2*exp(-0.5*sum((xxt+c).^2,1))/sqrt(2*pi)^op.dim;
gTRUEFig=-(p1*(xxt-c).*exp(-0.5*sum((xxt-c).^2,1))/sqrt(2*pi)^op.dim...
    +p2*(xxt+c).*exp(-0.5*sum((xxt+c).^2,1))/sqrt(2*pi)^op.dim)./pdfFig;

figure(1); subplot(122);
plot(xxt,gTRUEFig,'-r','LineWidth',1.5); 
title('A Mixture of two Gaussians','FontSize',15); 
axis([-5 5 -3 3]); legend('LSLDG','True'); 
axis square; set(gca,'FontSize',15); hold off;