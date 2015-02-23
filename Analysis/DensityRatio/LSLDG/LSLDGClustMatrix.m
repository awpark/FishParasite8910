function Y=LSLDGClustMatrix(X,op,C)
%
% Clustering via mode seeking based on the LSLDG estimator with a matrix-form update.
%

narginchk(2,3);

if nargin==2; [~,theta,C,sigma]=LSLDG(X,op);
else [~,theta,C,sigma]=LSLDG(X,op,C);
end

Y=zeros(size(X));

iter=0; criteria=1;
while criteria
    iter=iter+1;
        
    if iter==1; Yold=X; end;
    
    Gau=permute(exp(-sum(bsxfun(@minus,permute(Yold,[1,3,2]),C).^2,1)...
        /(2*sigma^2)),[2,3,1]);
    Y=((theta.*C)*Gau)./(theta*Gau); 
    
    if op.disp && iter<=20
        figure(1); subplot(122); scatter(Y(1,:),Y(2,:),'ko'); 
        axis([-5 5 -5 5]); axis square; 
        title(['Iter=' num2str(iter)]); set(gca,'FontSize',15); 
    end

    % convergence criteria 
    d=sqrt(sum((Y-Yold).^2,1)); d=max(d(d<inf),[],2);
    criteria = (d > op.tol) & (iter < op.maxiter);

    Yold=Y;        
end

% For irregular data points
infind=find(isinf(sum(abs(Y),1))==1);
if ~isempty(infind); Y(:,infind)=NaN; end;

% Assigning clusters to irregular data points
if sum(isnan(sum(Y,1)))~=0; Y=assignCluster(Y,X); end;
