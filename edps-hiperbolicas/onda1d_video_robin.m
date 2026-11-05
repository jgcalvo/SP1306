% Ecuacion de onda en 1D: video de la condicion de Robin en x=0
%
%   d_nu u + lambda u_t = 0,    r = (1 - lambda c)/(1 + lambda c)
%
% Tres filas: lambda*c = 1/3 (r = +1/2), lambda*c = 1 (r = 0, frontera absorbente)
% y lambda*c = 3 (r = -1/2). Los casos limite son Neumann (lambda = 0, r = +1) y
% Dirichlet (lambda -> infinito, r = -1).
%
% Escribe onda_robin.mp4 en el directorio actual.

AZ=[.10 .35 .62]; RJ=[.78 .25 .12]; GRIS=[.62 .62 .62];

c = 1;
g = @(x) (x>=1 & x<=2).*sin(pi*(x-1));      % mismo pulso que el otro video
x = linspace(0,6,1601);

% u(x,t) = 1/2 g(x+ct) + 1/2 g(x-ct) + (r/2) g(ct-x)
sol = @(r,t) .5*g(x+c*t) + .5*g(x-c*t) + (r/2)*g(c*t-x);

lamc = [1/3 1 3];                            % lambda*c
rr   = (1-lamc)./(1+lamc);
etiq = {sprintf('\\lambda c = 1/3   \\rightarrow   r = %+.1f   (refleja la mitad, sin invertir)',rr(1)), ...
        sprintf('\\lambda c = 1      \\rightarrow   r = %.0f      (frontera absorbente: la onda se va)',rr(2)), ...
        sprintf('\\lambda c = 3      \\rightarrow   r = %+.1f   (refleja la mitad, invertida)',rr(3))};
fprintf('r = %+.3f  %+.3f  %+.3f\n',rr);

fps=30; Tfin=5; nf=round(15*fps); tt=linspace(0,Tfin,nf);
nhold=round(0.8*fps);

f=figure('Visible','off','Color','w','Position',[0 0 1280 860]);
ax=gobjects(1,3);
ax(1)=axes('Position',[.075 .695 .895 .255]);
ax(2)=axes('Position',[.075 .400 .895 .255]);
ax(3)=axes('Position',[.075 .105 .895 .255]);

vw=VideoWriter('onda_robin','MPEG-4');
vw.FrameRate=fps; vw.Quality=95; open(vw);

for k=[ones(1,nhold) 1:nf]
  t=tt(k);
  for j=1:3
    u=sol(rr(j),t);
    cla(ax(j)); hold(ax(j),'on');
    plot(ax(j),[0 6],[0 0],'-','Color',[.82 .82 .82],'LineWidth',.6);
    plot(ax(j),[0 6],[.5 .5],':','Color',[.85 .85 .85],'LineWidth',.8);
    plot(ax(j),[0 6],[-.5 -.5],':','Color',[.85 .85 .85],'LineWidth',.8);
    plot(ax(j),x,u,'-','Color',AZ,'LineWidth',3.0);
    plot(ax(j),[0 0],[-1.15 1.15],'-','Color',RJ,'LineWidth',4);
    hold(ax(j),'off');
    set(ax(j),'XLim',[0 6],'YLim',[-1.2 1.2],'YTick',[-1 -.5 0 .5 1], ...
              'FontSize',13,'Box','off','TickDir','out','LineWidth',.8);
    text(ax(j),.18,-.92,etiq{j},'FontSize',15,'Color',RJ);
  end
  xlabel(ax(3),'x','FontSize',15);
  text(ax(1),5.8,.95,sprintf('t = %4.2f',t),'FontSize',17, ...
       'HorizontalAlignment','right','Color',[.25 .25 .25]);
  writeVideo(vw,getframe(f));
end
close(vw); close(f);
fprintf('listo: onda_robin.mp4  (%d cuadros, %.1f s, %d fps)\n', ...
        nf+nhold, (nf+nhold)/fps, fps);
