function [x,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D(x0,nit,par,varargin)
% construction a permuation matrix and then reconstruct the MATRIX such
%Every row is sample!!!
% that \|X0-PX\|+par1\|X\|_{tv,l1} 
% An illustrative example:
%x0=[1+0.1*rand(10,1);2+0.1*rand(10,1);0.1*rand(40,1)];
%[x,P,prim,Dual]=permu_TVL1_Secular(x0(truth)/mean(x0),500,[2,0.1]);
% 结果是把样本按列平均化了
%% Parameters
if nargin < 4
    tol = 1e-4;
end
    
      
beta=10;
gamma=(1+sqrt(5))/2*beta;
gamma=5*beta;
alpha = par(1);

[m,n ] = size(x0);
D =gallery('circul',[1,-1,zeros(1,m-2)]);
P =  eye(m);
t=0;
CF=zeros(nit,1);

x=x0;
y=D*x0;    %=nabla x0

lambda=0.5*ones(size(y));

%% Computation of the cost function at start

CF(1)=sum(x(:)-x0(:));
B1= inv(gamma/2*D'*D+eye(m))/2; 
norm_tv = tvnorm(x);
 
%% Core of the algorithm: ADM methods
for i=2:nit
    %Define y=(w,z)
    %L(x,y,lambda) = ||w||_1 + chi_K(z) + <lambda, Ax-y> + beta/2 ||Ax-y||^2_2
    %% First step : x^{k+1}=argmin_x (L(x,y^{k},lambda^{k}))
    %x = B1*(x0'*P-D'*lambda+gamma*D'*y+x0*P');
    %x = B1*(P'*x0+P*x0'-D'*(lambda-gamma*y));
       x = B1*(2*P'*x0-D'*(lambda-gamma*y));
    minx = min(x(:)); maxx = max(x(:));
    %x = (x-minx)/(maxx-minx);
     %% Third step: update P     %A'= B'*P';
    P = estimate_permuation(x0,x);
%% Second step : y^{k+1}=argmin_x (L(x^{k+1},y,lambda^{k}))
    yzprec = y;
%     xx1=(D*x/2+1./gamma.*x);
%     for i=1:size(y,2)
%         nxx=abs(xx1(:,i));nxx(nxx==0)=1;
%     beta1= gamma/2/alpha;
%      y(:,i)=xx1(:,i)-min(1./beta1,nxx).*xx1(:,i)./nxx;
%     end
   E = D*x/2+1./gamma.*x;
   y = sign(E) .* max(abs(E)-gamma/2/alpha, 0);

     %% Third step: update P     %A'= B'*P';
    P = estimate_permuation(x0,x);
    %% Third step : lambda^{k+1}=lambda^k+gamma (Nabla x^{k+1}-y{k+1})
    lambda= lambda+beta*(D*x-y);
    %% Fourth step : adjusting gamma

   
  DualRes(i) =gamma*(norm(yzprec-y,2)); % primal residual
  PrimRes(i) = norm( D*x-y,2); % dual residual
  norm_tv(i) =tvnorm(x);
 fprintf('Iter %d: PrimRes = %f, DualRes = %f, gamma = %f \n',i,PrimRes(i),DualRes(i),gamma);
 fprintf('TV Norm = %f \n',log(norm_tv(i)));  
   %th = norm(x(:))*tol;
  %% Fourth step : adjusting beta
 
 if  abs((PrimRes(i)-PrimRes(i-1)))/PrimRes(i) < tol && abs((DualRes(i)-DualRes(i-1)))/DualRes(i) < tol
        break
    end
     if PrimRes(i)>0.5*DualRes(i)
       % gamma = 1.2*gamma; %beta=gamma;
        beta = 1.2*beta;
    elseif DualRes(i)>0.5*PrimRes(i)
       % gamma = gamma/1.2; %beta=gamma;
        beta = beta/1.2;
     else
     end  
     
      if ~mod(i,100)
       % subplot(nit/100,1,i/100); imagesc(x);   axis off;       title(strcat('Iter = ',num2str(i)));
        tempx{i/100}= x;
      end
end
%     if norm_tv(i)>norm_tv(i-1)
%         gamma = .5*gamma;
%     end
%     x = x/max(x(:));
%     y = y/max(y(:));
  
%smoothed= post_proc(x,x0);
end
function value =tvnorm(im)
d=zeros(size(im));
d(1:end-1,:)=im(2:end,:)-im(1:end-1,:);
d(end,:)=im(1,:)-im(end,:);
value = sum(abs(d(:)));
end

 


