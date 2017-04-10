function sx = TCP_postprocessing(x)
[C,indCluster]=max(x,[],1);
 [indClusterSorted,ind]=sort(indCluster);
sx=x(:,ind);