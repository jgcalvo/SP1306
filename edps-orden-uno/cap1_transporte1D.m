% Ecuacion de transporte, 1D
b = .1;                % velocidad
u = @(x,t) sin(x-b*t); % solucion exacta
% 2d plot
xx = linspace(-10,10);
tt = linspace(0,10);
[x,t] = meshgrid(xx,tt);
surf(x,t,u(x,t))
xlabel('x')
ylabel('t')
% 1d plot
figure
axis([-10 10 -1 1])
for j = 1:5:numel(tt)    % se grafican algunas curvas
    plot(xx,u(xx,tt(j)))
    title(num2str(tt(j)))
    pause(.1)            % para observar el cambio en tiempo
end