

for i =2359:2369

    p(:,:,i)=p(:,:,i)./repmat(sum(p(:,:,i)),8,1); % make probabilities add up to 1 across classes for each exemplar

    for j=1:max(classes)
        avg(:,j)=mean(p(:, 1+(j-1)*5:j*5, i),2);
        avg_shift(:,j)=circshift(avg(:,j),4-j); % align to position 4 as correct 
    end;

    figure; hold on;
    plot(mean(avg_shift,2),'b');
    plot(mean(avg_shift,2)+std(avg_shift,0,2)./sqrt(40),'b--');
    plot(mean(avg_shift,2)-std(avg_shift,0,2)./sqrt(40),'b--');

    plot([1 8],[1/8 1/8],'r'); % chance;
    set(gcf,'color','w');
    xlabel('correct angle - classified angle');
    ylabel('probability');
    %set(gca,'YLim',[0 0.2]);

end