function permu_mat = estimate_permuation(A,B)
%A = PB each row is sample
%A'= B'*P';
 sample_size = size(A,1);
% A = rand(sample_size,25);
% truth= randperm(sample_size);
% B(1:sample_size,:) =A(truth,:) ;
 D = pdist2(A,B,'minkowski');  
% D = pdist2(A,B,'cityblock');
[Assign,cost] = munkres(D);
permu_mat = zeros(sample_size);
for i=1:sample_size
    permu_mat(i,Assign(i))=1;
end

