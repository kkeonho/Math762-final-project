function [ ] = movie(x,N,Lx,Ly)

% [n,m] = size(x);


for i = 1:N
%     if n == 1
%         plot(x(i),0,'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
%         axis([-1 1 -1 1])
%         grid on
%         drawnow
%     else
        figure(5)
        scatter(x(:,2*i-1),x(:,2*i),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
%         plot(x(i,:),0,'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
%         axis([-Lx*0.7 Lx*0.7 -Ly*0.55 Ly*0.55])
%         axis([-Lx/4 5*Lx/4 -Ly*0.1 Ly*1.1])
        grid on
        drawnow
%     end
%     pause(0.001);
end
end

