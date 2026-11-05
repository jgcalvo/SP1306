% Ecuacion de onda en 1D: video de la reflexion en x=0
%
% Compara Dirichlet (extremo fijo) con Neumann (extremo libre). En gris punteado
% van las dos mitades viajeras de la formula de d'Alembert; su suma es la curva
% azul. Asi se ve que en Neumann las dos se superponen en la pared y dan el doble
% de la amplitud incidente, y que en Dirichlet se cancelan.
%
% Escribe onda_reflexion.mp4 en el directorio actual.
% La figura se mantiene oculta a proposito: getframe captura mejor asi.

AZ=[.10 .35 .62]; RJ=[.78 .25 .12]; GRIS=[.62 .62 .62];

c   = 1;
g   = @(x) (x>=1 & x<=2).*sin(pi*(x-1));     % pulso en [1,2]
go  = @(x) sign(x).*g(abs(x));               % extension impar  (Dirichlet)
ge  = @(x) g(abs(x));                        % extension par    (Neumann)
x   = linspace(0,6,1601);

fps=30; Tfin=5; nf=round(15*fps); tt=linspace(0,Tfin,nf);
nhold=round(0.8*fps);                        % pausa inicial

f=figure('Visible','off','Color','w','Position',[0 0 1280 720]);
ax1=axes('Position',[.075 .555 .895 .345]);
ax2=axes('Position',[.075 .105 .895 .345]);

vw=VideoWriter('onda_reflexion','MPEG-4');
vw.FrameRate=fps; vw.Quality=95; open(vw);

for k=[ones(1,nhold) 1:nf]
  t=tt(k);
  pD=.5*go(x-c*t); mD=.5*go(x+c*t);   uD=pD+mD;
  pN=.5*ge(x-c*t); mN=.5*ge(x+c*t);   uN=pN+mN;

  dibuja(ax1,x,uD,pD,mD,AZ,RJ,GRIS,'Dirichlet:  extremo fijo,  u(0,t)=0');
  dibuja(ax2,x,uN,pN,mN,AZ,RJ,GRIS,'Neumann:  extremo libre,  u_x(0,t)=0');
  xlabel(ax2,'x','FontSize',15);
  text(ax1,5.75,.95,sprintf('t = %4.2f',t),'FontSize',17,'FontName','Helvetica',...
       'HorizontalAlignment','right','Color',[.25 .25 .25]);

  writeVideo(vw,getframe(f));
end
close(vw); close(f);

fprintf('listo: onda_reflexion.mp4  (%d cuadros, %.1f s, %d fps)\n', ...
        nf+nhold, (nf+nhold)/fps, fps);

function dibuja(ax,x,u,p,m,AZ,RJ,GRIS,etiq)
  cla(ax); hold(ax,'on');
  plot(ax,x,p,'--','Color',GRIS,'LineWidth',1.1);     % componente ->
  plot(ax,x,m,'--','Color',GRIS,'LineWidth',1.1);     % componente <-
  plot(ax,x,u,'-','Color',AZ,'LineWidth',3.0);        % solucion
  plot(ax,[0 6],[0 0],'-','Color',[.8 .8 .8],'LineWidth',.6);
  plot(ax,[0 0],[-1.15 1.15],'-','Color',RJ,'LineWidth',4);
  hold(ax,'off');
  set(ax,'XLim',[0 6],'YLim',[-1.2 1.2],'YTick',[-1 -.5 0 .5 1], ...
         'FontSize',13,'Box','off','TickDir','out','LineWidth',.8);
  text(ax,.18,-.95,etiq,'FontSize',15,'FontName','Helvetica','Color',RJ);
end
