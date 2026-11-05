% Ecuacion de onda en 1D: video de la cuerda pulsada
%
% Dato inicial triangular g(x) = min(x, 1-x), velocidad inicial nula, un periodo
% completo t en [0, 2L/c]. A la izquierda el perfil exacto, calculado con la
% extension impar 2-periodica; a la derecha los modos n = 1, 3, 5 con su amplitud
% real, que hace visible el decaimiento A_n ~ n^(-2).
%
% Escribe cuerda_pulsada.mp4 en el directorio actual.

AZ=[.10 .35 .62]; RJ=[.78 .25 .12]; GRIS=[.62 .62 .62];
COL=[AZ; .85 .55 .15; .35 .60 .35];           % modos n=1,3,5

L=1; c=1;
gtri = @(x) min(x,1-x).*(x>=0 & x<=1);
x    = linspace(0,1,1201);
ns   = [1 3 5];
An   = (4./(ns.^2*pi^2)).*sin(ns*pi/2);
fprintf('A_1=%.4f  A_3=%.4f  A_5=%.4f\n',An);

fps=30; Tfin=2; nf=round(13*fps); tt=linspace(0,Tfin,nf);
nhold=round(0.8*fps);

f=figure('Visible','off','Color','w','Position',[0 0 1280 620]);
ax1=axes('Position',[.055 .145 .415 .74]);
ax2=axes('Position',[.555 .145 .415 .74]);

vw=VideoWriter('cuerda_pulsada','MPEG-4');
vw.FrameRate=fps; vw.Quality=95; open(vw);

for k=[ones(1,nhold) 1:nf]
  t=tt(k);
  u = .5*(per2odd(x-c*t,gtri)+per2odd(x+c*t,gtri));

  % --- panel izquierdo: la cuerda ---------------------------------------
  cla(ax1); hold(ax1,'on');
  plot(ax1,x,gtri(x),'-','Color',[.85 .85 .85],'LineWidth',1.2);
  plot(ax1,[0 1],[0 0],'-','Color',[.8 .8 .8],'LineWidth',.6);
  plot(ax1,x,u,'-','Color',AZ,'LineWidth',3.2);
  plot(ax1,[0 0],[-.6 .6],'-','Color',RJ,'LineWidth',4);
  plot(ax1,[1 1],[-.6 .6],'-','Color',RJ,'LineWidth',4);
  hold(ax1,'off');
  set(ax1,'XLim',[-.012 1.012],'YLim',[-.62 .62],'YTick',-.5:.25:.5, ...
          'FontSize',13,'Box','off','TickDir','out','LineWidth',.8);
  xlabel(ax1,'x','FontSize',15);
  title(ax1,'La cuerda','FontSize',16,'FontWeight','normal');
  text(ax1,.97,.54,sprintf('t = %4.2f',t),'FontSize',16,...
       'HorizontalAlignment','right','Color',[.25 .25 .25]);
  if max(abs(u))<.02
    text(ax1,.5,-.44,'cuerda plana: toda la energía es cinética', ...
         'FontSize',14,'HorizontalAlignment','center','Color',RJ);
  end

  % --- panel derecho: los modos -----------------------------------------
  cla(ax2); hold(ax2,'on');
  plot(ax2,[0 1],[0 0],'-','Color',[.8 .8 .8],'LineWidth',.6);
  for j=1:3
    n=ns(j);
    plot(ax2,x,An(j)*sin(n*pi*x)*cos(n*pi*c*t/L),'-', ...
         'Color',COL(j,:),'LineWidth',2.4);
  end
  plot(ax2,[0 0],[-.6 .6],'-','Color',RJ,'LineWidth',4);
  plot(ax2,[1 1],[-.6 .6],'-','Color',RJ,'LineWidth',4);
  hold(ax2,'off');
  set(ax2,'XLim',[-.012 1.012],'YLim',[-.62 .62],'YTick',-.5:.25:.5, ...
          'FontSize',13,'Box','off','TickDir','out','LineWidth',.8);
  xlabel(ax2,'x','FontSize',15);
  title(ax2,'Los tres primeros modos (amplitud real)','FontSize',16,'FontWeight','normal');
  lg=legend(ax2,{'','n = 1','n = 3','n = 5'},'FontSize',13,'Box','off', ...
            'Location','southeast');
  lg.ItemTokenSize=[18 10];

  writeVideo(vw,getframe(f));
end
close(vw); close(f);
fprintf('listo: cuerda_pulsada.mp4  (%d cuadros, %.1f s, %d fps)\n', ...
        nf+nhold, (nf+nhold)/fps, fps);

function y=per2odd(x,g)
  z=mod(x,2); y=zeros(size(z));
  m1=z<=1;  y(m1)=g(z(m1));
  m2=~m1;   y(m2)=-g(2-z(m2));
end
