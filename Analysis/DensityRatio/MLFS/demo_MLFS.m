clear all;

rand('state',0);
randn('state',0);

dataset=1;
switch dataset
  case 1 % classification
   n=200;
   X1=[randn(2,n/2).*repmat([1;2],[1 n/2])+repmat([-3;0],[1 n/2])];
   X2=[randn(2,n/2).*repmat([1;2],[1 n/2])+repmat([ 3;0],[1 n/2])];
   X=[X1 X2];
   Y=[-ones(1,n/2) ones(1,n/2)];
   y_type=1; % classification
   colormap_type=prism;
  case 2 % regression
   n=1000;
   X=(rand(2,n)*2-1)*10;
   Y=sin(X(1,:)/10*pi);
   y_type=0; % regression
   colormap_type=hsv;
end

[feature_order,MIh]=MLFS(X,Y,y_type);

feature_order
MIh

%%%%%%%%%%%%%%%%%%%%%% Displaying original 2D data
figure(dataset)
clf
hold on

colormap(colormap_type)
set(gca,'FontName','Helvetica')
set(gca,'FontSize',12)
tmp=(feature_order==1);
scatter3(X(1,:),X(2,:),Y,100,Y,'filled');
h=plot([-tmp(1) tmp(1)]*100,[-tmp(2) tmp(2)]*100,'k-','LineWidth',4);
axis equal
axis([-10 10 -10 10])
title('Original 2D data and important feature found by MLFS')

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 12 12]);
print('-depsc',sprintf('MLFS%g',dataset))
  
