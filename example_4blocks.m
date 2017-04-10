function example_4blocks

m=100;
A= eye(m)+0.05*rand(m);
A = toeplitz([1  0 1 0 1 0 1 0  zeros(1,84) 1  0 1 0 1 0 1 0]);
%A = eye(m);
w= [A ;
    0.01*rand(150,m)];;
A = rand(m); [Q,R] = qr(A);

w = [Q-min(Q(:));
    0.01*rand(150,m)];; 
 k=100;
  l = 0:25:m;
 h = 0.1*rand(m);
 for i=1:length(l)-1
     h(l(i)+1:l(i+1),l(i)+1:l(i+1)) = i+0.05*rand(25);
 end
 %h = [h; 0.01*randn(40,size(h,2))]; %h = [h, randn(size(h,1),40)];
 h1 = h;
order2=randperm(size(h1,1));
h2 = h1(order2,:);
OP = zeros(m);
for i=1:m
    OP(i,order2(i))=1;
end
% h2 = OP*h1;

%optimal_norm = norm(OP*h1-h2)+ 500*tvnorm(h1);
% alpha = 1e+2;
% sigma  =0.5;
% nit = 600;
% noise = sigma*randn(size(h2));
% [x,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(h2+noise,nit,alpha);
% figure; subplot(3,1,1); imagesc(x);
% subplot(3,1,2); imagesc(h2+noise);
% subplot(3,1,3); imagesc(h1+noise);

% Experiment 1 to demonstrate the robustness of noisy feature

P = test_noisy_feature(h2,1000);
P = test_noisy_feature(h2,100);
P = test_noisy_feature(h2,10);

% Experiment 2 to demonstrate the  
draw_figure(0.01,h2);
draw_figure(0.1,h2);
draw_figure(0.5,h2)
draw_figure(1,h2)
draw_figure(2,h2)
 
end

function P = test_noisy_feature(h2,ndim)
%ndim = 500;
sigma = [0.02 0.2 2];
alpha = 1e+2;
nit = 300;
threshold = 0.99;
h3 = [h2 0.1*randn(size(h2,1),ndim)];
rp_h3 = h3(:,randperm(size(h3,2)));
figure;
for i=1:length(sigma)
nh3 = rp_h3+sigma(i)*randn(size(h3));
%Dx = FD_pca(nh3,threshold);
[x,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(nh3,nit,alpha);
residue = nh3-P*x;
snr_value = mean((nh3(:)).^2)/mean(residue(:).^2);
subplot(3,3,i); imagesc(nh3); axis off; title(strcat('SNR = ',num2str(snr_value)));
 subplot(3,3,i+3); imagesc(P'*h3); axis off; %title('Before Removing Noisy Features');
subplot(3,3,i+6); imagesc(P'*h2); axis off; %title('After Removing Noisy Features'); 
end
print(gcf,'-depsc', strcat('Exp2_Dim', num2str(ndim),'.eps'));
 

end

function value =tvnorm(im)
d=zeros(size(im));
d(1:end-1,:)=im(2:end,:)-im(1:end-1,:);
d(end,:)=im(1,:)-im(end,:);
value = sum(abs(d(:)));
end

function Dx = FD_pca(ndata,threshold)
[pc,score,latent,tsquare] = princomp(ndata);
dimension = find(cumsum(latent)./sum(latent)>threshold);
tranMatrix = pc(:,1:dimension(1));
Dx = ndata*tranMatrix;
 
end

function draw_figure(sigma,h2)
 alpha = 1e+2;
nit = 500;
noise = sigma*randn(size(h2));
[x,P,PrimRes,norm_tv,tempx]=permu_TVL1_Secular_2D_v2(h2+noise,nit,alpha);
residue = h2+noise-P*x;
snr_value = mean((h2(:)+noise(:)).^2)/mean(residue(:).^2);
figure; subplot(4,1,1); imagesc(h2+noise); axis off; title(strcat('SNR = ',num2str(snr_value)));
subplot(4,1,2); imagesc(tempx{1});  axis off;  title('Iter = 100');
subplot(4,1,3); imagesc(tempx{3});  axis off;   title('Iter = 300');
subplot(4,1,4); imagesc(x);   axis off;   title('Final');
print(gcf,'-depsc', strcat('Fig4simulated_SNR',num2str(snr_value),'.eps'));
end