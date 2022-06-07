close all
clear;

% define function
% f   = @(x) x.^2 + sin(pi*x);
% df  = @(x) 2*x + pi*cos(pi*x);
% ddf = @(x) 2 - pi*pi*sin(pi*x);
f   = @(x) 1/4 * x.^2 + sin(pi*x);
df  = @(x) 1/2 * x + pi*cos(pi*x);
ddf = @(x) 1/2 - pi*pi*sin(pi*x);

% define plot domain/space
X   = -3*pi:0.001:3*pi;
Y   = f(X);
% dX  = df(X);
% ddX = ddf(X);

% define problem hyperparams
te = 0.0001;
tc = 0.00001;
lri = 1;
momentum = 0.9;
x0 = 3*pi;

figure('Position', [0 0 512 384]);
file = 'anim.gif';
s = 5;

x = x0;
y = f(x);
Xs = repmat(x, 1, s);
Ys = repmat(y, 1, s);
dx = 0;

stop = 0;
i = 0;
turn = 0;
max_turn = 3;
while (~stop)
    % solve the problem
    i = i + 1;
    lr = lri / 1;
    dx  = dx * momentum + df(x) * (1 - momentum);
    ddx = ddf(x);
    x = x - lr * (dx);
    y = f(x);
    
    % show up
    Xs = [Xs(2:s) x];
    Ys = [Ys(2:s) y];
    clf('reset');
    plot(X, Y, 'Color', 'b');
    hold on;
%     plot(X, ddX, 'Color', 'm');
    for j = 1:s
        pl = scatter(Xs(j), Ys(j), 'o', 'MarkerEdgeColor', 'r');
        pl.MarkerEdgeAlpha = j^2/s^2;
    end
    grid on;
    set(gcf,'color','w');
    xlabel(['lr = ' num2str(lr) '; momentum = ' num2str(momentum) '; x_0 = ' num2str(x0) '; dd = ' num2str(ddx)]);
    frame(i) = getframe(gcf);
%     pause(0.1);
    
    % stopping condition
    stop = (sum(abs(Xs(s) - Xs)) < s * tc);
    if (stop)
        turn = turn + 1;
        
    end
    if (turn < max_turn)
        stop = 0;
    end
end
hold off;

% save animation to file
for i = 1:length(frame)
    [imind,cm] = rgb2ind(frame(i).cdata,256);
    if i == 1 
        imwrite(imind,cm,file,'gif', 'Loopcount',inf,'DelayTime',0.0333); 
    else 
        imwrite(imind,cm,file,'gif','WriteMode','append','DelayTime',0.0333); 
    end 
end