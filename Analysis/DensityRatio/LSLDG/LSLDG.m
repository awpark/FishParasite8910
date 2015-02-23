function [g,theta,C,sigma,lambda,cind]=LSLDG(X,op,C)
%
% Estimating the Gradient of a log-density
% 
% X: (dim) by (sample) matrix
% op: options
% C: center points of Gaussian functions
%

narginchk(2,3);

if nargin==2
    cind=randperm(op.samples,op.bnum);
    C=X(:,cind);
else
    cind=[];
end

% Difference to centers (size: bnum by samles by dim)
XC_diff=repmat(permute(C,[2,3,1]),[1,op.samples,1])...
    -repmat(permute(X,[3,2,1]),[op.bnum,1,1]);

% Distance from centers (size: bnum by samples)
XC_dist=sum(XC_diff.^2,3);

% cross validation
cv_fold=(1:op.cvfold);
cv_split=floor((0:op.samples-1)*op.cvfold./op.samples)+1;
cv_index=cv_split(randperm(op.samples));

score_cv=zeros(length(op.sigma_list),length(op.lambda_list),length(cv_fold));
for sigma_index=1:length(op.sigma_list)
    sigma=op.sigma_list(sigma_index);   

    GauKer3D=repmat(exp(-XC_dist/(2*sigma^2)),[1,1,op.dim]);
    
    for k=cv_fold      
        psi_train=XC_diff(:,cv_index~=k,:).*GauKer3D(:,cv_index~=k,:);      
        phi_train=(1-XC_diff(:,cv_index~=k,:).^2/sigma^2).*GauKer3D(:,cv_index~=k,:);
      
        [G_train,H_train]=computeGH(psi_train,phi_train,op);
      
        psi_test=XC_diff(:,cv_index==k,:).*GauKer3D(:,cv_index==k,:);      
        phi_test=(1-XC_diff(:,cv_index==k,:).^2/sigma^2).*GauKer3D(:,cv_index==k,:);

        [G_test,H_test]=computeGH(psi_test,phi_test,op);

        for lambda_index=1:length(op.lambda_list)
            lambda=op.lambda_list(lambda_index);

            term1=0; term2=0;
            for dd=1:op.dim
                Gtmp_train=G_train(:,(dd-1)*op.bnum+1:dd*op.bnum);
                htmp_train=H_train(:,dd);
                
                Gtmp_test=G_test(:,(dd-1)*op.bnum+1:dd*op.bnum);
                htmp_test=H_test(:,dd);

                thetah=mldivide(Gtmp_train+lambda*eye(size(Gtmp_train)),htmp_train);
                    
                term1=term1+thetah'*Gtmp_test*thetah;
                term2=term2+thetah'*htmp_test;                    
            end % d                          
            score_cv(sigma_index,lambda_index,k)=term1-2*term2;
        end % lambda
    end % k
end % sigma

[score_cv_tmp,lambda_index]=min(mean(score_cv,3),[],2);
[~,sigma_index]=min(score_cv_tmp);
lambda=op.lambda_list(lambda_index(sigma_index));
sigma=op.sigma_list(sigma_index);
fprintf('sigma=%g, lambda=%g\n',sigma,lambda);

GauKer3D=repmat(exp(-XC_dist/(2*sigma^2)),[1,1,op.dim]);

psi=XC_diff.*GauKer3D;
phi=(1-XC_diff.^2/sigma^2).*GauKer3D;

[G,H]=computeGH(psi,phi,op);

theta=zeros(op.bnum,op.dim);
for dd=1:op.dim
    Gtmp=G(:,(dd-1)*op.bnum+1:dd*op.bnum);
    htmp=H(:,dd);

    theta(:,dd)=mldivide(Gtmp+lambda*eye(size(Gtmp)),htmp);        
end

g=permute(sum(bsxfun(@times,theta,permute(psi,[1,3,2])),1),[2,3,1]);

theta=theta';

function [G,H]=computeGH(psi,phi,op)
%
% Computing G and H in obj func.
% 

narginchk(3,3);

G=zeros(op.bnum,op.bnum*op.dim);
for dd=1:op.dim
    G(:,(dd-1)*op.bnum+1:dd*op.bnum)...
        =psi(:,:,dd)*psi(:,:,dd)'/size(psi,2);
end

H=permute(mean(phi,2),[1,3,2]);
