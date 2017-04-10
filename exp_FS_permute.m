function exp_FS_permute 
addpath(genpath('D:\FangCloudSync\Matlab\FSLib_v4.0_2016'));
sigma=0.1;
[X,h2,h1,class]= create_data(sigma,1000);
nit = 500;
alpha = 1e+2;
dim = 100;
 
[ranking_cfs,ranking_llcfs,fea_MCFS, ranking_Lap]=exp(X,dim);
[x_cfs,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(X(:,ranking_Lap),nit,alpha);
 sx_cfs = TCP_postprocessing(x_cfs);
%figure; imagesc(x);
[x_llcfs,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(X(:,ranking_llcfs),nit,alpha);
 sx_llcfs = TCP_postprocessing(x_llcfs');
%figure; imagesc(x);
[x_MCFS,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(X(:,fea_MCFS),nit,alpha);
 sx_MCFS = TCP_postprocessing(x_MCFS');
%figure; imagesc(x);
[x_Lap,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(X(:,ranking_Lap),nit,alpha);
 sx_Lap = TCP_postprocessing(x_Lap');
%figure; imagesc(x);
figure; 
subplot(5,1,1); imagesc(X);   axis off; title('True');
subplot(5,1,2); imagesc(sx_llcfs);  axis off; %title(strcat('Raw (SNR = ',num2str(snr_value),')'));
subplot(5,1,3); imagesc(sx_cfs);   axis off; title('Sparse NMF');
subplot(5,1,4); imagesc(sx_MCFS);   axis off;  title('Quick Sort');
subplot(5,1,5); imagesc(sx_Lap); axis off; title('TCP');
%print(gcf,'-depsc', strcat('FD4methods_snr',num2str(snr_value),'.eps'));

end

function [rp_h3,h3,h1,class]= create_data(sigma,ndim)
addpath(genpath('D:\FangCloudSync\Matlab\work\nmfv1_4'))
m=100; 
  l = 0:25:m;
 va =[ 1 2 3 4 5]
 h = 0.1*rand(m); 
  class = ones(1,m);
 for i=1:length(l)-1
     h(l(i)+1:l(i+1),l(i)+1:l(i+1)) = va(i)+0.05*rand(25);
      class(l(i)+1:l(i+1)) = i;
 end
 
 %h = [h; 0.01*randn(40,size(h,2))]; %h = [h, randn(size(h,1),40)];
 h1 = h;
order2=randperm(size(h1,1));
class = class(order2);
h2 = h1(order2,:);
OP = zeros(m);
for i=1:m
    OP(i,order2(i))=1;
end
noise = sigma*randn(size(h2));
nh2 = h2+noise;
h3 = [nh2 0.1*randn(size(h2,1),ndim)];
rp_h3 = h3(:,randperm(size(h3,2)));

end

function [ranking_cfs,ranking_llcfs,fea_MCFS, ranking_Lap]=exp(X,dim)
%[ranking_cfs,ranking_llcfs,fea_MCFS,
%index_LaplacianScore]=exp_FS_permute(Dx,dim)£»

temp = cfs(X); ranking_cfs = temp(1:dim);
temp = llcfs(X); ranking_llcfs = temp(1:dim);
options = [];
[FeaIndex,FeaNumCandi] = MCFS_p(X, dim,options);
for i = 1:length(FeaNumCandi)
    fea_MCFS = FeaIndex{i};
   end
options = [];
options.Metric = 'Cosine';
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'Cosine';
W = constructW(X,options);
L = LaplacianScore(X,W);
[junk, temp] = sort(-L);
ranking_Lap =  temp(1:dim);

end