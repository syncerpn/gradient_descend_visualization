close all
clear;

% define function
% f   = @(x) x.^2 + sin(pi*x);
% df  = @(x) 2*x + pi*cos(pi*x);
% ddf = @(x) 2 - pi*pi*sin(pi*x);
f     = @(x, y) 1/4 * x.^2 + sin(pi*x) + y.^2 + sin(pi*y);
df_x  = @(x, y) 1/2 * x + pi*cos(pi*x);
ddf_x = @(x, y) 1/2 - pi*pi*sin(pi*x);
df_y  = @(x, y) 2 * y + pi*cos(pi*y);
ddf_y = @(x, y) 2 - pi*pi*sin(pi*y);


% define plot domain/space
X   = -1*pi:0.1:1*pi;
Y   = -1*pi:0.1:1*pi;
[Xm, Ym] = meshgrid(X, Y);
Zm = f(Xm, Ym);

% dX  = df(X);
% ddX = ddf(X);

% define problem hyperparams
tc = 0.0001;
lri = 0.1;
momentum = 0.98;
x0 = -0.5*pi;
y0 = 0.5*pi;

figure('Position', [0 0 512 384]);
file = 'anim.gif';
s = 5;

x = x0;
y = y0;
z = f(x, y);
Xs = repmat(x, 1, s);
Ys = repmat(y, 1, s);
Zs = repmat(z, 1, s);
dx = 0;
dy = 0;

stop = 0;
i = 0;
while (~stop)
    % solve the problem
    i = i + 1;
    lr = lri / 1;
    dx  = dx * momentum + df_x(x, y) * (1 - momentum);
    dy  = dy * momentum + df_y(x, y) * (1 - momentum);
    ddx = ddf_x(x,y);
    ddy = ddf_y(x,y);
    x = x - lr * (dx);
    y = y - lr * (dy);
    z = f(x, y);
    
    % show up
    Xs = [Xs(2:s) x];
    Ys = [Ys(2:s) y];
    Zs = [Zs(2:s) z];
    clf('reset');
    contour(Xm, Ym, Zm, 'LevelStep', 1);
    hold on;
    
    for j = 1:s
        pl = scatter(Xs(j), Ys(j), 'o', 'MarkerEdgeColor', 'r');
        pl.MarkerEdgeAlpha = j^2/s^2;
    end
    grid on;
%     xlabel(['lr = ' num2str(lr) '; momentum = ' num2str(momentum) '; y_0 = ' num2str(x0) '; dd = ' num2str(ddy)]);
    set(gcf,'color','w');
    frame(i) = getframe(gcf);
    
    % stopping condition
    stop = (sum(abs(Xs(s) - Xs)) < s * tc) && (sum(abs(Ys(s) - Ys)) < s * tc);
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