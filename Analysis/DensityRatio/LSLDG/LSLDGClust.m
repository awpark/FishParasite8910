function Y=LSLDGClust(X,op,C)
%
% Clustering via mode seeling based on the LSLDG estimator with a vector-form update.
%

narginchk(2,3);

if nargin==2; [~,theta,C,sigma]=LSLDG(X,op);
else [~,theta,C,sigma]=LSLDG(X,op,C);
end

Y=zeros(size(X));

for ii=1:op.samples
 
    Z=X(:,ii);
    
    iter=0; criteria=1; NaNcall=0;
    while criteria
        iter=iter+1;
    
        Gau=exp(-sum(bsxfun(@minus,Z,C).^2)'/(2*sigma^2));        
        
        Zold=Z;
        Z=((theta.*C)*Gau)./(theta*Gau); 
        
        % For an iregular data point
        if sum(isnan(sum(Z,1)))~=0 || sum(isinf(sum(Z,1)))~=0
            NaNcall=1; break;
        end
    
        % convergence criteria 
        d=sqrt(sum((Z-Zold).^2));
        criteria = (d > op.tol) & (iter < op.maxiter);
    end

    if NaNcall; Y(:,ii)=NaN; else Y(:,ii)=Z; end;
end
  
% Assigning clusters to irregular points
if sum(isnan(sum(Y,1)))~=0; Y=assignCluster(Y,X); end;
