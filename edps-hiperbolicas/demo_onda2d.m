%=========================================================================
%  DEMO: ecuacion de onda en 2D con elemento finito P1
%
%     u_tt = Delta u   en (0,1)^2,    u = 0 en la frontera,
%     u(x,0) = g,      u_t(x,0) = 0.
%
%  Con c = 1 y g(x) = sen(3 pi x1) sen(4 pi x2) la solucion exacta es
%
%     u(x,t) = sen(3 pi x1) sen(4 pi x2) cos(5 pi t),
%
%  porque Delta u = -(9+16) pi^2 u = -25 pi^2 u = u_tt.  Periodo: 2/5.
%
%  El guion hace tres cosas:
%     (1) arma M y A, y calcula la CFL del esquema explicito;
%     (2) corre las tres discretizaciones en tiempo y compara con la exacta;
%     (3) barre dt para ver exactamente donde explota cada una.
%  Es autocontenido: no necesita ningun otro archivo.
%=========================================================================
clear; clc
ANIMAR = true;          % poner false para saltarse la animacion

AZ=[.10 .35 .62]; RJ=[.78 .25 .12]; VD=[.25 .50 .30]; NA=[.85 .55 .15];

%% ---------- (1) malla, matrices y CFL -----------------------------------
h = 1/32;
[M,A,xy,node,elem,in] = mallaP1(h);
uex = @(t) sin(3*pi*xy(:,1)).*sin(4*pi*xy(:,2))*cos(5*pi*t);
U0  = uex(0);   V0 = zeros(size(U0));

lmax  = eigs(A,M,1,'largestabs');
dtCFL = 2/sqrt(lmax);
fprintf('h = 1/%d,  %d nodos interiores\n',round(1/h),numel(U0));
fprintf('lambda_max = %.4e      lambda_max*h^2 = %.2f\n',lmax,lmax*h^2);
fprintf('dt_CFL = 2/sqrt(lambda_max) = %.4e\n\n',dtCFL);

%% ---------- (2) las tres discretizaciones -------------------------------
%  Newmark-beta con (beta,gamma):
%     (0, 1/2)        = diferencias centradas  (explicito, condicional)
%     (1/4, 1/2)      = Crank-Nicolson         (implicito, incondicional)
%     (0.3025, 0.6)   = Newmark disipativo     (implicito, incondicional, orden 1)
T = 0.8;                                   % dos periodos
casos = { 'Centradas,  dt = 0.9 CFL'      , 0     , .5 , 0.9*dtCFL , AZ
          'Crank-Nicolson, dt = 0.9 CFL'  , .25   , .5 , 0.9*dtCFL , VD
          'Crank-Nicolson, dt = 10 CFL'   , .25   , .5 , 10 *dtCFL , NA
          'Newmark (.3025,.6), dt = 10 CFL', .3025, .6 , 10 *dtCFL , RJ };

fprintf('(2) evolucion hasta T = %.2f\n',T);
fprintf('    %-32s %10s %7s %12s %12s\n','metodo','dt','pasos','err final','E(T)/E(0)');
res = cell(size(casos,1),1);
for j = 1:size(casos,1)
    dt = casos{j,4};  nt = round(T/dt);  dt = T/nt;
    [Uh,Vh,tv] = newmark(M,A,U0,V0,dt,nt,casos{j,2},casos{j,3});
    err = zeros(1,nt+1);  E = zeros(1,nt+1);
    for k = 1:nt+1
        err(k) = norm(Uh(:,k)-uex(tv(k)),inf);
        E(k)   = .5*(Vh(:,k)'*M*Vh(:,k)) + .5*(Uh(:,k)'*A*Uh(:,k));
    end
    res{j} = struct('t',tv,'err',err,'E',E/E(1),'U',Uh);
    fprintf('    %-32s %10.3e %7d %12.3e %12.4f\n', ...
            casos{j,1},dt,nt,err(end),E(end)/E(1));
end

figure('Name','Error y energia','Position',[80 400 980 340]);
subplot(1,2,1)
for j=1:4, plot(res{j}.t,res{j}.err,'-','Color',casos{j,5},'LineWidth',1.7); hold on, end
grid on; xlabel('t'); ylabel('|| u_h(t) - u(t) ||_\infty'); xlim([0 T])
legend(casos(:,1),'Location','northwest','Box','off','FontSize',8)
title('Error contra la solucion exacta')
subplot(1,2,2)
for j=1:4, plot(res{j}.t,res{j}.E,'-','Color',casos{j,5},'LineWidth',1.7); hold on, end
yline(1,':','Color',[.6 .6 .6]); grid on; xlabel('t'); ylabel('E^n / E^0')
xlim([0 T]); ylim([0 1.15]); title('Energia discreta')

%% ---------- (3) estabilidad: barrido de dt ------------------------------
fprintf('\n(3) barrido de dt (todas hasta T = %.2f)\n',T);
fprintf('    %8s %12s %14s %14s %14s\n','dt/CFL','dt', ...
        'centradas','C-N','Newmark disip.');
for f = [0.5 0.9 0.99 1.0 1.01 1.05 1.2 5 50]
    dt = f*dtCFL;  nt = max(2,round(T/dt));  dt = T/nt;
    m = zeros(1,3);  pr = [0 .5; .25 .5; .3025 .6];
    for i = 1:3
        U = newmark(M,A,U0,V0,dt,nt,pr(i,1),pr(i,2));
        m(i) = norm(U(:,end),inf);
    end
    fprintf('    %8.2f %12.3e %14.3e %14.3e %14.3e\n',f,dt,m);
end
% umbral fino para el esquema explicito, en dos horizontes distintos
fprintf('\n    umbral empirico del explicito (criterio ||U||_inf > 10):\n');
for Tu = [T 5*T]
    lo = 0.90;  hi = 1.10;
    while hi-lo > 1e-4
        f = (lo+hi)/2;  dt = f*dtCFL;  nt = round(Tu/dt);
        U = newmark(M,A,U0,V0,Tu/nt,nt,0,.5);
        if norm(U(:,end),inf) > 10, hi = f; else, lo = f; end
    end
    fprintf('      mirando hasta T = %4.1f :  dt = %.4f dt_CFL\n',Tu,lo);
end
fprintf('    El umbral exacto es 1. Con horizonte corto sale un poco mayor\n');
fprintf('    porque justo encima de la CFL el modo inestable crece despacio\n');
fprintf('    y no le da tiempo; al alargar T converge a 1.\n');
fprintf('\n    Lecturas:\n');
fprintf('    - El explicito pasa de O(1) a desbordarse justo al cruzar la CFL.\n');
fprintf('    - Los implicitos nunca se desbordan, pero eso no los vuelve\n');
fprintf('      correctos: con dt = 50 CFL la amplitud es 0.35 y 0.03 en vez de\n');
fprintf('      ~1. C-N se desfasa y Newmark ademas amortigua.\n');
fprintf('    - Ojo con las dos ultimas filas: con dt tan grande el explicito da\n');
fprintf('      solo 2 pasos y no le alcanza para desbordarse, pero el numero\n');
fprintf('      que imprime tampoco significa nada.\n');

%% ---------- animacion ----------------------------------------------------
if ANIMAR
    j = 1;  Uh = res{j}.U;  tv = res{j}.t;    % centradas a 0.9 CFL
    figure('Name','Evolucion','Position',[80 60 900 380]);
    for k = 1:4:numel(tv)
        subplot(1,2,1); pinta(node,elem,in,Uh(:,k));
        title(sprintf('u_h,   t = %.3f',tv(k)))
        subplot(1,2,2); pinta(node,elem,in,uex(tv(k)));
        title('exacta')
        drawnow
    end
end

%=========================================================================
%  funciones locales
%=========================================================================
function [M,A,xy,node,elem,in] = mallaP1(h)
%MALLAP1  Matrices de masa y rigidez P1 en (0,1)^2, solo nodos interiores.
    x = 0:h:1;  n = numel(x);
    [X,Y] = ndgrid(x,x);  node = [X(:) Y(:)];
    id = @(i,j) (j-1)*n + i;
    elem = zeros(2*(n-1)^2,3);  k = 0;
    for j = 1:n-1
        for i = 1:n-1
            k=k+1; elem(k,:) = [id(i,j)   id(i+1,j)   id(i+1,j+1)];
            k=k+1; elem(k,:) = [id(i,j)   id(i+1,j+1) id(i,j+1)  ];
        end
    end
    NT = size(elem,1);  N = size(node,1);
    ve = zeros(NT,2,3);                       % ve(:,:,i): lado opuesto al vertice i
    ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
    ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
    ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
    area = .5*(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
    ii=[]; jj=[]; sA=[]; sM=[];
    for a = 1:3
        for b = 1:3
            ii = [ii; elem(:,a)];  jj = [jj; elem(:,b)];
            sA = [sA; dot(ve(:,:,a),ve(:,:,b),2)./(4*area)];
            sM = [sM; area*((a==b)+1)/12];
        end
    end
    A = sparse(ii,jj,sA,N,N);   M = sparse(ii,jj,sM,N,N);
    bd = node(:,1)==0 | node(:,1)==1 | node(:,2)==0 | node(:,2)==1;
    in = ~bd;  A = A(in,in);  M = M(in,in);  xy = node(in,:);
end

function [Uh,Vh,tv] = newmark(M,A,U0,V0,dt,nt,be,ga)
%NEWMARK  Familia de Newmark-beta para M U'' + A U = 0, guardando la historia.
    U = U0;  V = V0;  Ac = M\(-A*U);
    K = M + be*dt^2*A;  R = chol(K);
    Uh = zeros(numel(U0),nt+1);  Vh = Uh;  Uh(:,1) = U0;  Vh(:,1) = V0;
    tv = (0:nt)*dt;
    for k = 1:nt
        Ut = U + dt*V + (dt^2/2)*(1-2*be)*Ac;      % predictor
        Un = R\(R'\(M*Ut));
        if be > 0, An = (Un-Ut)/(be*dt^2); else, An = M\(-A*Un); end
        V  = V + dt*((1-ga)*Ac + ga*An);           % corrector
        U  = Un;  Ac = An;
        Uh(:,k+1) = U;  Vh(:,k+1) = V;
        if ~all(isfinite(U)), Uh(:,k+1:end) = Inf; Vh(:,k+1:end) = Inf; return, end
    end
end

function pinta(node,elem,in,Uint)
%PINTA  Dibuja una funcion P1 dada por sus valores en los nodos interiores.
    u = zeros(size(node,1),1);  u(in) = Uint;
    trisurf(elem,node(:,1),node(:,2),u,'EdgeColor','none','FaceColor','interp');
    axis([0 1 0 1 -1.05 1.05]); clim([-1 1]); view(-37.5,30); shading interp
    xlabel('x_1'); ylabel('x_2')
end
