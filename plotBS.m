function plotBS(X,Y,r,K)
% Plot BS Deployment

for i = 1:length(X);
    xTemp = X(i,1);
    yTemp = Y(i,1);
    %plot(xTemp,yTemp,'r.');              % Plot only BS points
    %circle(xTemp,yTemp,r);               % Plot circular cell borders
    [xunit,yunit] = hexagon(xTemp,yTemp,r);
    plot(xunit, yunit,'k');               % Plot hexagonal cell borders 
    colors = {'r','b','g','y','m','c',[150 100 190]./255};
    fill(xunit,yunit,colors{mod(i,K)+1}); 
    text(xTemp,yTemp,[' ' num2str(i) ' '],'Color','k','Margin',0.1,...
        'FontSize',10,'HorizontalAlignment','center');    % Plot BS numbers
    axis equal, axis([0.1,0.9,0.1,0.9]), zoom on
    hold on
end
end