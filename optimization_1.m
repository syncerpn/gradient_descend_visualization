% close all
clf('reset');

clear;

% define function
% f   = @(x) x.^2 + sin(pi*x);
% df  = @(x) 2*x + pi*cos(pi*x);
% ddf = @(x) 2 - pi*pi*sin(pi*x);
f   = @(x) 1/4 * x.^2 + sin(pi*x);
df  = @(x) 1/2 * x + pi*cos(pi*x);
ddf = @(x) 1/2 - pi*pi*sin(pi*x);

% define plot domain/space
X   = -9:0.001:9;
Y   = f(X);
dX  = df(X);
ddX = ddf(X);

% define problem hyperparams
tc = 0.0001;
lri = 1;
momentum = 0.9;
x0 = 9;

% figure('Position', [0 0 512 384]);
file = 'anim.gif';
s = 2;

x = x0;
y = f(x);
dx = df(x);
ddx = ddf(x);
Xs = repmat(x, 1, s);
Ys = repmat(y, 1, s);
dXs = repmat(dx, 1, s);
ddXs = repmat(ddx, 1, s);
dxu = 0;

stop = 0;
i = 0;
while (~stop)
    % solve the problem
    i = i + 1;
    lr = lri / sqrt(i);
    dxu  = dxu * momentum + dx * (1 - momentum);
    x = x - lr * (dx);
    y = f(x);
    dx = df(x);
    ddx = ddf(x);
    
    % show up
    Xs = [Xs(2:s) x];
    Ys = [Ys(2:s) y];
    dXs = [dXs(2:s) dx];
    ddXs = [ddXs(2:s) ddx];
    
    clf('reset');
    plot(X, Y, 'Color', 'b');
    hold on;
    grid on;
    plot(X, dX, 'Color', 'c');
    plot(X, ddX, 'Color', 'm');
    for j = 1:s
        pl = scatter(Xs(j), Ys(j), 'o', 'MarkerEdgeColor', 'r');
        pl.MarkerEdgeAlpha = j^2/s^2;
        
        pld = scatter(Xs(j), dXs(j), 'x', 'MarkerEdgeColor', 'r');
        pld.MarkerEdgeAlpha = j^2/s^2;
        
        pldd = scatter(Xs(j), ddXs(j), 'x', 'MarkerEdgeColor', 'r');
        pldd.MarkerEdgeAlpha = j^2/s^2;
        
        pll = plot([Xs(j), Xs(j)], [dXs(j), ddXs(j)], '-', 'Color', 'r');
    end
    set(gcf,'color','w');
    xlabel(['lr = ' num2str(lr) '; momentum = ' num2str(momentum) '; x_0 = ' num2str(x0) '; dd = ' num2str(ddx)]);
%     frame(i) = getframe(gcf);
    pause(0.1);
    
    % stopping condition
    stop = (sum(abs(Xs(s) - Xs)) < s * tc);
end
hold off;

% save animation to file
% for i = 1:length(frame)
%     [imind,cm] = rgb2ind(frame(i).cdata,256);
%     if i == 1 
%         imwrite(imind,cm,file,'gif', 'Loopcount',inf,'DelayTime',0.0333); 
%     else 
%         imwrite(imind,cm,file,'gif','WriteMode','append','DelayTime',0.0333); 
%     end 
% end