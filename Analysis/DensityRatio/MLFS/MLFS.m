function [feature_order,MIh]=MLFS(x,y,y_type,feature_group,sigma_list,b,fold)
%
% Maximum Likelihood Feature Selection
%
% Order input features according to their dependency on outputs
%
% Usage:
%       [feature_order,MIh]=MLFS(x,y,y_type,feature_group,sigma_list,b)
%
% Input:
%    x             : dx by n input sample matrix
%    y             : dy by n output sample matrix
%    y_type     : if y_type=1, delta kernel is used for y; 
%                 otherwise (or empty) Gaussian kernel is used.
%    feature_group : dx-dimensional vector indicating group indices;
%                    if dx features are grouped, their group indices are specified.
%                    Default is [1:dx].
%    sigma_list    : (candidates of) Gaussian width
%                    If sigma_list is a vector, one of them is selected by cross validation.
%                    If sigma_list is a scalar, this value is used without cross validation
%                    If sigma_list is empty/undefined, Gaussian width is chosen from
%                    some default canditate list by cross validation
%    b             : number of Gaussian centers (if empty, b=200 is used);
%
% Output:
%    feature_order : Feature orders 
%                    feature_order(1)  : most dependent
%                    feature_order(end): least dependent
%    MIh           : Estimated mutual information for each feature
%
% (c) Taiji Suzuki, Department of Mathematical Informatics, The University of Tokyo, Japan. 
%     Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     s-taiji@stat.t.u-tokyo.ac.jp
%     sugi@cs.titech.ac.jp,

[dx n] =size(x);
[dy ny] =size(y);
if n~=ny
    error('number of samples of x and y must be the same!!!')
end

if nargin<3 || isempty(y_type)
  y_type=0;
end

if nargin < 4 || isempty(feature_group)
    feature_group = [1:dx]; 
end
if length(feature_group) ~= dx
    error('the length of feature_group must be equal to the dimension of x!!!')    
end

if nargin < 5  || isempty(sigma_list)
  sigma_list=logspace(-2,2,9);
end

if nargin < 6 || isempty(b)
  b=200;
end

if nargin < 7 || isempty(fold)
  fold=5;
end

    
global mytol
mytol = 10^(-15);

group_list=unique(feature_group);

for group_index=1:length(group_list)
  z=x(find(feature_group==group_list(group_index)),:);

  %Gaussian centers are randomly chosen from samples
  rand_index=randperm(n);
  u=z(:,rand_index(1:b));
  v=y(:,rand_index(1:b));
  
  dy=size(y,1);
  Phi_tmp =GaussBasis_sub([z;y],[u;v])';
  Phiz_tmp=GaussBasis_sub(z,u)';
  if y_type==0
    Phiy_tmp=GaussBasis_sub(y,v)';
  end
  
  if length(sigma_list)==1 
    sigma_chosen=sigma_list; %no cross validation
    score_cv=-inf;
  else %Choose sigma by cross validation
    fold_index=[1:fold];
    cv_index=randperm(n);
    cv_split=floor([0:n-1]*fold./n)+1;
    scores_cv=zeros(length(sigma_list),1);
    
    for sigma_index=1:length(sigma_list)
      sigma=sigma_list(sigma_index);
      Phiz_sigma=GaussBasis(Phiz_tmp,sigma);
      if y_type
        Phiy_sigma=DeltaBasis(y,v);
      else
        Phiy_sigma=GaussBasis(Phiy_tmp,sigma);
      end
      Phi_sigma=Phiz_sigma.*Phiy_sigma;
      
      for i=1:fold
        cv_index_tmp=cv_index(cv_split==i);
        bbz(:,i)=sum(Phiz_sigma(:,cv_index_tmp),2);
        bby(:,i)=sum(Phiy_sigma(:,cv_index_tmp),2);
        bbzy(:,i)=sum(Phiz_sigma(:,cv_index_tmp).*Phiy_sigma(:,cv_index_tmp),2);
        n_cv(i)=length(cv_index_tmp);
      end
      
      for i=1:fold
        cv_index_tr=fold_index(fold_index~=i);
        n_cv_tr=sum(n_cv(cv_index_tr));
        bb=(sum(bbz(:,cv_index_tr),2).*sum(bby(:,cv_index_tr),2)...
            -sum(bbzy(:,cv_index_tr),2))/((n_cv_tr-1)*n_cv_tr);
        alphah_cv=KLIEP_learning(bb,Phi_sigma(:,cv_index(cv_split~=i))');
        wh_cv=alphah_cv'*Phi_sigma(:,cv_index(cv_split==i));
        scores_cv(sigma_index)=scores_cv(sigma_index)+mean(log(wh_cv(wh_cv>mytol)))/fold;
      end %fold
    end %sigma_index
    
    scores_cv(find(isnan(scores_cv)))=-Inf;
    [score_cv,sigma_index_cv]=max(scores_cv);
    sigma_chosen=sigma_list(sigma_index_cv(1));
  end %length(sigma_list)==1 
  
  %%%%%%%%%%%%%%%% Computing the final solution `MIh'
  Phiz=GaussBasis(Phiz_tmp,sigma_chosen);
  if y_type
    Phiy=DeltaBasis(y,v);
  else
    Phiy=GaussBasis(Phiy_tmp,sigma_chosen);
  end
  Phi=Phiz.*Phiy;
  bb=(sum(Phiz,2).*sum(Phiy,2) - sum(Phiz.*Phiy,2))/(n*(n-1));
  alphah=KLIEP_learning(bb,Phi');
  wh=alphah'*Phi;
  MIh(group_index)=mean(log(wh(wh>mytol)));
end %group_index

[dummy Ind]=sort(MIh,'descend');
feature_order=Ind;
