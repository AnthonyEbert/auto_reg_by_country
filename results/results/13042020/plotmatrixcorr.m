function [S,ax] = plotmatrixcorr(data,labels,titlestr,CI,subset,highlight)

if nargin > 3
data = data(:, subset);
labels = labels(subset);
CI.c95 = CI.c95(:,subset);


I = highlight;
else
    I = [];
end

N = length(labels);

[S,ax,bx] = plotmatrix(data);
[rho,pval] = corr(data);
for i=1:N-1
    
    if nargin > 3
        hold on
        for j=i+1:N
            %cla(ax(j,i));
            %axes();
            errorbar(ax(j,i),data(:,i),data(:,j),CI.c95(:,j),CI.c95(:,j),...
                CI.c95(:,i),CI.c95(:,i),'o','Marker','none','CapSize',0);
            if (length(I) > 0)
                errorbar(ax(j,i),data(I,i),data(I,j),CI.c95(I,j),CI.c95(I,j),...
                CI.c95(I,i),CI.c95(I,i),'or','CapSize',0);
                xlim([min(data(:,i)),max(data(:,i))]);
                ylim([min(data(:,j)),max(data(:,j))]);
            end
            if (i > 1)
                set(ax(j,i),'YTickLabel',[]);
            end
            if (j < N)
                set(ax(j,i),'XTickLabel',[]);
            end
        end
        hold off
    end
    for j=i+1:N
        cla(ax(i,j));
        axes(ax(i,j));
        %text(0,0,['$\rho=',num2str(rho(i,j),'(\text{p-val }=',pval(i,j),')$')]);
        axis([0,1,0,1]);
        str = ['$',sprintf('%0.2f',rho(i,j)),'$'];
        if (rho(i,j) < -0.004)
            text(0.5,0.5,str,'interpreter','latex','color','red');
        elseif (rho(i,j) > 0.004)
            str = ['$+',sprintf('%0.2f',rho(i,j)),'$'];
            text(0.5,0.5,str,'interpreter','latex','color','blue');
        else
            rho(i,j) = 0;
            text(0.5,0.5,str,'interpreter','latex','color','black');
        end
            %['$\rho=',num2str(rho(i,j)),'(\text{p-val }=',pval(i,j),')$']
    end
end
for i=1:N
    set(ax(i,1).YLabel,'string',['$',labels{i},'$']);
    set(ax(N,i).XLabel,'string',['$',labels{i},'$']);
end
axes(bx)
title(titlestr)