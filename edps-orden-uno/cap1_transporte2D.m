% Ecuacion de transporte, 2D
b = [10 2]; % velocidad
u = @(x,y,t) sin(sqrt((x-b(1)*t).^2+(y-b(2)*t).^2)); % solucion
% 2d plot
xx = linspace(-20,20,1e3);
tt = linspace(0,15,150);
[x,y] = meshgrid(xx,xx);
% guardar video 
writerObj = VideoWriter('transporte.avi');
writerObj.FrameRate = 20;
open(writerObj);
% grafique para algunos valores de t
figure, axis([-20 20 -20 20 -1 1])
xlabel('x')
ylabel('t')
for j = 1:2:numel(tt)
    surf(x,y,u(x,y,tt(j))), shading interp
    title(num2str(tt(j)))
    axis([-20 20 -20 20 -1 1])
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    pause(.1)
end
close(writerObj);