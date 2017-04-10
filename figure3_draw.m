function heatmap_zcs()

load nbReal1Finnal
load results4SRC3_city_a_1000g_5

r=randperm(size(data_M,1));   %生成行数的随机数序列
rd=data_M(r, :);    %将100个cells的RD信号随机排列


%RD信号热图, 论文中Figure4(a)
H1 = HeatMap(round(rd)-2,'Colormap',redbluecmap(5),'DisplayRange',2)
addYLabel(H1,'Cell Number')
addXLabel(H1,'Read Depth');
set(H1,'ColumnLabelsLocation','top')
p1 = plot(H1);
c1 = colorbar(p1);
set(c1,'TickLabels',[ 0 1 2 3 4])
set(c1,'YTick',[-2 -1 0 1 2 3 ])
print(gcf,'-depsc2',  'Fig4RD6.eps')
%TCP结果, 论文中Figure4(b)
H2 = HeatMap(round(x)-2,'Colormap',redbluecmap(5),'DisplayRange',2)
% addYLabel(H2,'Cell Number')
% addXLabel(H2,'Read Depth');
% p2 = plot(H2);
% c2 = colorbar(p2);
% set(c1,'YTick',[-2 -1 0 1 2 3 ])
%print(gcf,'-depsc',  'Fig4TCP.eps')

%对100个cells的RD信号按行hierarchical clustering, 论文中Figure4(c)


clu1=clustergram(round(data_M)-2,'Cluster',1,'RowPDist','cityblock',...
            'RowLabels',PloidyInfo,'colormap',redbluecmap(5),'DisplayRange',4);
   
end

