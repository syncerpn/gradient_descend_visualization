close all
clear;

% define function
f   = @(x) x.^2;
df  = @(x) 2*x;

% define plot domain/space
X   = -4:0.01:4;
Y   = f(X);

% define problem hyperparams
tc  = 0.0001;
lri = 0.1;
momentum = 0.9;
x0 = 3.9;

figure('Position', [1000 600 512 384]);
file = 'anim_0.gif';
s = 5;

x = x0;
y = f(x);
Xs = repmat(x, 1, s);
Ys = repmat(y, 1, s);
dx = 0;

stop = 0;
i = 0;
y_best = inf;

while (~stop)
    % solve the problem
    i = i + 1;
    lr = lri;
    dx  = dx * momentum + df(x) * (1 - momentum);
    x = x - lr * (dx);
    y = f(x);
    
    y_best = min(y, y_best);
    
    % show up
    Xs = [Xs(2:s) x];
    Ys = [Ys(2:s) y];
    clf('reset');
    plot(X, Y, 'Color', 'b');
    hold on;
    
    for j = 1:s
        pl = scatter(Xs(j), Ys(j), 'o', 'MarkerEdgeColor', 'r');
        pl.MarkerEdgeAlpha = j^2/s^2;
    end
    grid on;
    set(gcf,'color','w');
    xlabel(['lr = ' num2str(lr) '; momentum = ' num2str(momentum) '; x_0 = ' ...
        num2str(x0) '; y = ' num2str(y) '; y_{best} = ' num2str(y_best)]);
    frame(i) = getframe(gcf);
    pause(0.0001);
    
    % stopping condition
    stop = (sum(abs(Xs(s) - Xs)) < s * tc);
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