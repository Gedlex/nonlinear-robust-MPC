function handle = plot_tube(z,M,delta,col,lw)
% plot ellipse in form ||x-z||_M = delta

phi = linspace(0,2*pi,100);
x_c = [cos(phi); sin(phi)];
x = chol(M)\x_c*delta;

handle = plot(x(1,:) + z(1),x(2,:) + z(2),'Color',col,'LineWidth',lw);
end