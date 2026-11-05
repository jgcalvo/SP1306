% Ecuacion de onda en 1D: figuras estaticas
%
% Genera tres archivos PDF en el directorio actual:
%   onda_reflexion.pdf     reflexion de un pulso en x=0, Dirichlet contra Neumann
%   onda_periodicidad.pdf  pulso rebotando en [0,1]: la solucion es 2L/c periodica
%   cuerda_pulsada.pdf     sumas parciales de la serie y perfil de la cuerda
%
% Imprime ademas las verificaciones numericas: u(0,t)=0 con Dirichlet, la
% duplicacion de amplitud con Neumann, la periodicidad exacta, los coeficientes
% A_n contra cuadratura, y la cuerda plana en t = L/(2c).

set(0,'DefaultAxesFontSize',9);
AZ=[.10 .35 .62]; RJ=[.78 .25 .12]; GR=[.45 .45 .45];

%% ---------- Figura 1: reflexion en la semirrecta  (Dirichlet vs Neumann)
g  = @(x) (x>=1 & x<=2).*sin(pi*(x-1));          % pulso positivo en [1,2]
go = @(x) sign(x).*g(abs(x));                     % extension impar
ge = @(x) g(abs(x));                              % extension par
x  = linspace(0,6,3001);  ts=[0 1.5 2.5 3.5];  c=1;

f=figure('Visible','on','Position',[0 0 980 340]);
for k=1:4
  t=ts(k);
  uD=.5*(go(x-c*t)+go(x+c*t));   uN=.5*(ge(x-c*t)+ge(x+c*t));
  subplot(2,4,k)
    plot(x,uD,'-','Color',AZ,'LineWidth',1.6); hold on
    plot([0 0],[-1.1 1.1],'-','Color',RJ,'LineWidth',2);
    yline(0,'-','Color',GR,'LineWidth',.4);
    axis([0 6 -1.15 1.15]); box off
    title(sprintf('t = %.1f',t),'FontWeight','normal')
    if k==1, ylabel({'Dirichlet','u(0,t)=0'},'Color',RJ); end
  subplot(2,4,4+k)
    plot(x,uN,'-','Color',AZ,'LineWidth',1.6); hold on
    plot([0 0],[-1.1 1.1],'-','Color',RJ,'LineWidth',2);
    yline(0,'-','Color',GR,'LineWidth',.4);
    axis([0 6 -1.15 1.15]); box off; xlabel('x')
    if k==1, ylabel({'Neumann','u_x(0,t)=0'},'Color',RJ); end
end
exportgraphics(gcf,'onda_reflexion.pdf','ContentType','vector');
fprintf('--- semirrecta ---\n');
fprintf('Dirichlet: u(0,t) max sobre t en [0,4] = %.2e\n', ...
   max(abs(arrayfun(@(t).5*(go(-t)+go(t)),linspace(0,4,401)))));
uN15=.5*(ge(x-1.5)+ge(x+1.5));
fprintf('Neumann  : max u(.,1.5) = %.4f   (amplitud incidente 0.5)\n',max(uN15));
uD35=.5*(go(x-3.5)+go(x+3.5));
fprintf('Dirichlet: min u(.,3.5) = %.4f   (pulso invertido)\n',min(uD35));

%% ---------- Figura 2: periodicidad en el intervalo [0,1]
L=1; c=1;
g2 = @(x) (x>=.1 & x<=.4).*sin(pi*(x-.1)/.3).^2;        % pulso suave
gt = @(x) per2odd(x,g2);                                 % impar 2-periodica
xx = linspace(0,1,2001); ts2=[0 .3 .6 1.0 1.6 2.0];
f=figure('Visible','on','Position',[0 0 980 330]);
for k=1:6
  t=ts2(k); u=.5*(gt(xx-c*t)+gt(xx+c*t));
  subplot(2,3,k)
    plot(xx,u,'-','Color',AZ,'LineWidth',1.6); hold on
    plot([0 0],[-1 1],'-','Color',RJ,'LineWidth',2);
    plot([1 1],[-1 1],'-','Color',RJ,'LineWidth',2);
    yline(0,'-','Color',GR,'LineWidth',.4);
    axis([0 1 -1.12 1.12]); box off
    title(sprintf('t = %.1f',t),'FontWeight','normal')
    if k>3, xlabel('x'); end
end
exportgraphics(gcf,'onda_periodicidad.pdf','ContentType','vector');
u0=.5*(gt(xx)+gt(xx)); u2=.5*(gt(xx-2)+gt(xx+2));
fprintf('--- intervalo ---\nperiodicidad: max|u(.,2L/c)-u(.,0)| = %.3e\n',max(abs(u2-u0)));

%% ---------- Figura 3: la cuerda pulsada
gtri=@(x) min(x,1-x).*(x>=0 & x<=1);
gtp = @(x) per2odd(x,gtri);
N=1:400; An=(4./(N.^2*pi^2)).*sin(N*pi/2);
Anq=zeros(size(N));
for n=N, Anq(n)=2*integral(@(s)gtri(s).*sin(n*pi*s),0,1,'AbsTol',1e-13); end
fprintf('--- cuerda pulsada ---\nmax|An formula - An cuadratura| = %.2e\n',max(abs(An-Anq)));

f=figure('Visible','on','Position',[0 0 940 300]);
subplot(1,2,1)
  xs=linspace(0,1,1201);
  plot(xs,gtri(xs),'-','Color',[0 0 0],'LineWidth',2.2); hold on
  est={'-','-','-'}; Ns=[1 3 11]; col=[RJ; .85 .55 .15; AZ];
  for j=1:3
    S=zeros(size(xs));
    for n=1:Ns(j), S=S+An(n)*sin(n*pi*xs); end
    plot(xs,S,est{j},'Color',col(j,:),'LineWidth',1.4);
  end
  legend({'g(x)','N=1','N=3','N=11'},'Location','south','FontSize',8,'Box','off');
  xlabel('x'); title('Sumas parciales en t=0','FontWeight','normal'); box off
subplot(1,2,2)
  tt=[0 .15 .3 .5]; sty={'-','-','-','--'};
  for j=1:4
    u=.5*(gtp(xs-tt(j))+gtp(xs+tt(j)));
    plot(xs,u,sty{j},'Color',AZ*(1-.22*(j-1))+.22*(j-1),'LineWidth',1.7); hold on
  end
  yline(0,'-','Color',GR,'LineWidth',.4);
  legend({'t=0','t=0.15','t=0.3','t=0.5'},'Location','northeast','FontSize',8,'Box','off');
  axis([0 1 -.05 .55]); xlabel('x'); box off
  title('Perfil de la cuerda','FontWeight','normal');
exportgraphics(gcf,'cuerda_pulsada.pdf','ContentType','vector');
u05=.5*(gtp(xs-.5)+gtp(xs+.5));
fprintf('en t=L/(2c)=0.5 la cuerda esta plana: max|u| = %.2e\n',max(abs(u05)));
% serie truncada contra reflexion
for Nt=[11 51 201]
  S=zeros(size(xs)); t=.3;
  for n=1:Nt, S=S+An(n)*sin(n*pi*xs)*cos(n*pi*t); end
  ur=.5*(gtp(xs-t)+gtp(xs+t));
  fprintf('N=%3d: max|serie - reflexion| en t=0.3 = %.3e\n',Nt,max(abs(S-ur)));
end

function y=per2odd(x,g)       % extension impar y 2-periodica de g:[0,1]->R
  z=mod(x,2); y=zeros(size(z));
  m1=z<=1;  y(m1)=g(z(m1));
  m2=~m1;   y(m2)=-g(2-z(m2));
end
