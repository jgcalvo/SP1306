%MEF_CALOR2D  Ecuacion del calor en 2D con elementos P1 y el metodo theta.
%
%   Resuelve
%
%       u_t = Lap(u)          en Omega x (0,T],  Omega = (0,1)x(0,1)
%       u = 0                 sobre la frontera
%       u(x,y,0) = u0(x,y)
%
%   con elementos finitos P1 sobre una triangulacion, avanzando en el tiempo
%   con Euler hacia adelante, Euler hacia atras y Crank-Nicolson.
%   Es la version en 2D de MEF_calor1D.m.
%
%   FORMULACION DEBIL.  Multiplicando por v que se anula en la frontera e
%   integrando por partes el termino de segundo orden,
%
%       (u_t, v) + (grad u, grad v) = (f, v)      para toda v
%
%   Al escribir u_h = sum_j u_j(t) phi_j aparecen las DOS matrices que devuelve
%   assembling.m: la de masa viene del termino temporal y la de rigidez del
%   termino espacial,
%
%       M u'(t) + A u(t) = b
%
%   En el problema estacionario M no se usa; aqui es la que carga la derivada
%   en el tiempo, y por eso el ensamblaje del capitulo sirve tal cual.
%
%   EL METODO THETA.  Los tres esquemas son el mismo, evaluando el termino A u
%   como una combinacion del tiempo viejo y el nuevo:
%
%       M (u^{n+1}-u^n)/dt + A (theta u^{n+1} + (1-theta) u^n) = b
%
%   es decir
%
%       (M + theta dt A) u^{n+1} = (M - (1-theta) dt A) u^n + dt b
%
%       theta = 0    Euler hacia adelante (explicito)
%       theta = 1    Euler hacia atras    (implicito)
%       theta = 1/2  Crank-Nicolson
%
%   El metodo "explicito" TAMBIEN pide resolver un sistema en cada paso, porque
%   M no es diagonal. Lo explicito es el tratamiento del operador espacial, no
%   el costo del paso. (La forma usual de volverlo genuinamente explicito es
%   condensar la masa; no se hace aqui, pero mas abajo se dice cuanto cambia la
%   estabilidad.)
%
%   En los tres casos la matriz de la izquierda NO cambia de un paso al otro.
%   Se factoriza UNA sola vez, fuera del ciclo, y cada paso es solo una
%   sustitucion mas un producto con la matriz de la derecha.
%
%   ESTABILIDAD.  Si lambda son los valores propios de A v = lambda M v, cada
%   modo se multiplica en cada paso por
%
%       g(z) = (1 - (1-theta) z) / (1 + theta z),      z = dt*lambda > 0
%
%   Pidiendo |g| <= 1 se obtiene la condicion
%
%       z (1 - 2 theta) <= 2
%
%   de modo que para theta >= 1/2 no hay restriccion ninguna: el implicito y
%   Crank-Nicolson son incondicionalmente estables. Para theta < 1/2 hay que
%   respetar dt <= 2/((1-2 theta) lambda_max), que con theta = 0 es el
%   dt <= 2/lambda_max de Euler hacia adelante. El paso lo fija el valor propio
%   MAS GRANDE, que corresponde al modo mas oscilatorio: justo el que la malla
%   resuelve peor y que no aporta nada fisico.
%
%   En 1D habia formula cerrada (lambda_max = 12/h^2, dt* = h^2/6). En 2D no:
%   sobre esta triangulacion la matriz de masa tiene estrella de 7 puntos y no
%   es un producto de Kronecker. Maximizando el simbolo se obtiene
%
%       lambda_max -> 25.856/h^2        de donde      dt* -> h^2/12.93
%
%   y con la masa condensada seria lambda_max -> 8/h^2, o sea dt* = h^2/4. La
%   masa consistente es unas 3.2 veces mas restrictiva, el mismo fenomeno que
%   en 1D (donde el factor era exactamente 3). Como no hay formula general para
%   mallas arbitrarias, dtCritico calcula lambda_max con eigs.
%
%   Del otro lado del espectro, lambda_min -> 2*pi^2, que es la tasa de
%   decaimiento exacta del dato inicial que se usa aqui.
%
%   UNA ADVERTENCIA SOBRE CRANK-NICOLSON.  Ser incondicionalmente estable no es
%   lo mismo que amortiguar. Cuando z -> infinito,
%
%       theta = 1   (implicito):       g -> 0      los modos rigidos se apagan
%       theta = 1/2 (Crank-Nicolson):  g -> -1     sobreviven, alternando signo
%
%   Esa es la diferencia entre A-estable y L-estable. Con datos iniciales
%   suaves, como el de aqui, no se nota; con datos rugosos o discontinuos,
%   Crank-Nicolson con dt grande produce oscilaciones espurias que no decaen.
%
%   ENSAMBLAJE.  La funcion local ensamblar es assembling.m: mismo formato de
%   tripletes, 3x3 = 9 entradas por triangulo, y la identidad
%
%       int_T grad(phi_i).grad(phi_j) = (v_i . v_j) / (4|T|)
%
%   donde v_i es el lado opuesto al vertice i. La matriz de masa sale del mismo
%   ciclo con area*((i==j)+1)/12.
%
%   ORDENES.  Orden 2 en h para los tres. En el tiempo, orden 1 para Euler
%   (hacia adelante o hacia atras) y orden 2 para Crank-Nicolson.

clear; clc

%% --- 1. Datos del problema -------------------------------------------
T   = 0.05;                            % tiempo final
u0f = @(x,y) sin(pi*x).*sin(pi*y);     % dato inicial
f   = @(x,y) zeros(size(x));           % sin fuente: u_t = Lap(u)
uex = @(x,y,t) exp(-2*pi^2*t).*sin(pi*x).*sin(pi*y);

% Con u0 = sin(pi x) sin(pi y) la solucion exacta es e^{-2 pi^2 t} u0: la forma
% se conserva y solo decae la amplitud. En 2D la tasa es 2*pi^2, el doble que
% en 1D, asi que conviene un T mas corto para que quede algo que medir.

n  = 32;                               % subdivisiones por lado
dt = T/200;                            % paso de tiempo

[node,elem] = squaremesh_basico(n);
ue = uex(node(:,1),node(:,2),T);

fprintf('n = %d  ->  %d nodos, %d triangulos\n', n, size(node,1), size(elem,1));
fprintf('dt = %.2e (%d pasos);  dt critico del explicito: %.3e\n', ...
    dt, round(T/dt), dtCritico(node,elem,0));
fprintf('\n   metodo            err_max       err_L2\n');
for s = {{'Euler atras',1.0},{'Crank-Nicolson',0.5}}
    nombre = s{1}{1}; th = s{1}{2};
    uk = resolverTheta(node,elem,dt,T,u0f,f,th);
    fprintf('   %-16s  %.4e   %.4e\n', nombre, norm(uk-ue,inf), ...
        errorL2(node,elem,uk,@(x,y) uex(x,y,T)));
end
% Cuidado al leer esa tabla: Euler hacia atras sale con MENOS error que
% Crank-Nicolson, y no es que sea mejor. Su error temporal y el espacial tienen
% signos opuestos -- el esquema decae de menos, la malla decae de mas -- y para
% este dt se cancelan de casualidad. El error de Crank-Nicolson, en cambio, ya
% es practicamente el error espacial puro, que es la precision real de la
% discretizacion. La cancelacion desaparece al cambiar dt o la malla; la
% seccion 4 lo muestra con claridad.

%% --- 2. Evolucion en el tiempo ----------------------------------------
instantes = [0 0.005 0.02 0.05];
figure
for k = 1:numel(instantes)
    uk = resolverCN(node,elem,dt,instantes(k),u0f,f);
    subplot(2,2,k)
    trisurf(elem, node(:,1), node(:,2), uk, 'EdgeColor','none');
    axis([0 1 0 1 0 1]); colorbar
    xlabel('x'); ylabel('y')
    title(sprintf('t = %g', instantes(k)))
end
sgtitle('Crank-Nicolson')

%% --- 3. Convergencia en el tiempo -------------------------------------
%  Para ver el orden temporal limpio NO se compara contra la solucion exacta,
%  sino contra una referencia calculada en la MISMA malla con un dt diminuto.
%  Asi el error espacial, que es identico en las dos, se cancela exactamente, y
%  lo que queda es puro error temporal. Es la unica forma de ver el orden 2 de
%  Crank-Nicolson: contra la solucion exacta su error cae por debajo del error
%  espacial de inmediato y la medicion se satura.
fprintf('\n===== ORDEN EN EL TIEMPO  (n = %d, referencia: CN con dt = T/40000) =====\n', n);
uref = resolverCN(node,elem,T/40000,T,u0f,f);

nsT = [10 20 40 80 160 320];
dtT = T./nsT;
errBEt = zeros(size(nsT));  errCNt = zeros(size(nsT));
fprintf('   pasos      dt       err BE    orden     err CN    orden\n');
for k = 1:numel(nsT)
    errBEt(k) = norm(resolverBE(node,elem,dtT(k),T,u0f,f) - uref, inf);
    errCNt(k) = norm(resolverCN(node,elem,dtT(k),T,u0f,f) - uref, inf);
    if k == 1
        o1 = NaN; o2 = NaN;
    else
        o1 = log2(errBEt(k-1)/errBEt(k)); o2 = log2(errCNt(k-1)/errCNt(k));
    end
    fprintf('%8d  %.2e  %.3e  %6.2f   %.3e  %6.2f\n', ...
        nsT(k), dtT(k), errBEt(k), o1, errCNt(k), o2);
end

% El explicito solo existe por debajo de dt*, muy a la izquierda de ese rango
dtc  = dtCritico(node,elem,0);
nsFE = [1000 2000 4000];
dtFE = T./nsFE;
errFEt = arrayfun(@(d) norm(resolverFE(node,elem,d,T,u0f,f)-uref,inf), dtFE);
fprintf('Euler explicito (solo donde es estable, dt* = %.2e):\n', dtc);
fprintf('%8d  %.2e  %.3e\n', [nsFE; dtFE; errFEt]);

% --- la prediccion teorica, exacta para este dato inicial ---------------
%  El dato inicial es un modo propio, asi que en la malla se propaga como un
%  solo modo discreto: la solucion semidiscreta es exp(-lambda_h t) v_h y el
%  metodo theta la reemplaza por g(dt*lambda_h)^n v_h. El error es entonces
%
%      |g(dt*lambda_h)^n - exp(-lambda_h T)|
%
%  sin constantes libres. Desarrollando para dt chico se recuperan los terminos
%  de libro, con lambda = lambda_h:
%
%      Euler atras:     e ~ exp(-lambda T) lambda^2 dt   T / 2
%      Crank-Nicolson:  e ~ exp(-lambda T) lambda^3 dt^2 T / 12
%
%  Ojo que ahi va lambda_h y no 2*pi^2: usar el valor continuo deja un sesgo
%  del orden de (lambda_h/2pi^2)^3 - 1, que en esta malla es medio por ciento.
lh   = lambdaMin(node,elem);
gBE  = @(z) 1./(1+z);
gCN  = @(z) (1-z/2)./(1+z/2);
gFE  = @(z) 1-z;
teoBEt = abs(gBE(dtT.*lh).^nsT  - exp(-lh*T));
teoCNt = abs(gCN(dtT.*lh).^nsT  - exp(-lh*T));
teoFEt = abs(gFE(dtFE.*lh).^nsFE - exp(-lh*T));
fprintf('lambda_h = %.6f  (2*pi^2 = %.6f)\n', lh, 2*pi^2);
fprintf('razon teoria/medicion:  BE en [%.4f, %.4f],  CN en [%.4f, %.4f]\n', ...
    min(teoBEt./errBEt), max(teoBEt./errBEt), ...
    min(teoCNt./errCNt), max(teoCNt./errCNt));

%% --- 4. Convergencia en el espacio ------------------------------------
%  Ahora contra la solucion exacta, con dt FIJO. Crank-Nicolson exhibe el orden
%  2 en h sin problema; Euler hacia atras no puede, porque su error temporal se
%  queda constante mientras el espacial baja, y los ordenes medidos se vuelven
%  absurdos (llegan a dar negativos cuando los dos errores, que tienen signos
%  opuestos, se cancelan). Ese contraste es la razon practica para preferir un
%  esquema de orden 2 en el tiempo.
fprintf('\n===== ORDEN EN EL ESPACIO  (dt = T/200 fijo) =====\n');
nsH = [4 8 16 32 64];
hH  = 1./nsH;
errBEh = zeros(size(nsH));  errCNh = zeros(size(nsH));
fprintf('    n       h       err BE    orden     err CN    orden\n');
for k = 1:numel(nsH)
    [nd,el] = squaremesh_basico(nsH(k));
    uek = uex(nd(:,1),nd(:,2),T);
    errBEh(k) = norm(resolverBE(nd,el,T/200,T,u0f,f) - uek, inf);
    errCNh(k) = norm(resolverCN(nd,el,T/200,T,u0f,f) - uek, inf);
    if k == 1
        o1 = NaN; o2 = NaN;
    else
        o1 = log2(errBEh(k-1)/errBEh(k)); o2 = log2(errCNh(k-1)/errCNh(k));
    end
    fprintf('%5d  %.4f  %.3e  %6.2f   %.3e  %6.2f\n', ...
        nsH(k), hH(k), errBEh(k), o1, errCNh(k), o2);
end

%% --- 5. Graficos del error --------------------------------------------
figure

subplot(1,2,1)
% los simbolos son lo medido; las lineas, la prediccion teorica
loglog(dtT,  teoBEt, '-',  'Color',[0 0.447 0.741], 'LineWidth',1.2, ...
    'DisplayName','teoria (modo unico)')
hold on
loglog(dtT,  teoCNt, '-',  'Color',[0.851 0.325 0.098], 'LineWidth',1.2, ...
    'HandleVisibility','off')
loglog(dtFE, teoFEt, '-',  'Color',[0.929 0.694 0.125], 'LineWidth',1.2, ...
    'HandleVisibility','off')
loglog(dtT,  errBEt, 'o', 'MarkerSize',7, 'LineWidth',1.3, ...
    'Color',[0 0.447 0.741],     'DisplayName','Euler atras')
loglog(dtT,  errCNt, 's', 'MarkerSize',7, 'LineWidth',1.3, ...
    'Color',[0.851 0.325 0.098], 'DisplayName','Crank-Nicolson')
loglog(dtFE, errFEt, '^', 'MarkerSize',7, 'LineWidth',1.3, ...
    'Color',[0.929 0.694 0.125], 'DisplayName','Euler adelante')
xline(dtc, 'r-', 'dt^*', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off')
hold off
xlabel('dt'); ylabel('error temporal')
title(sprintf('Error vs dt  (n = %d):  pendientes 1 y 2', n))
legend('Location','southeast'); grid on

subplot(1,2,2)
loglog(hH, errBEh, 'o-', 'LineWidth',1.4, 'DisplayName','Euler atras')
hold on
loglog(hH, errCNh, 's-', 'LineWidth',1.4, 'DisplayName','Crank-Nicolson')
loglog(hH, 0.35*errCNh(1)*(hH/hH(1)).^2, 'k:', 'DisplayName','pendiente 2')
hold off
xlabel('h'); ylabel('error total')
title('Error vs h  (dt = T/200 fijo)')
legend('Location','southeast'); grid on

%% --- 6. La restriccion de estabilidad del explicito -------------------
%  Misma malla y mismo dt para los tres. Los dos incondicionalmente estables
%  funcionan en toda la columna; el explicito solo por debajo de dt*, y apenas
%  se pasa no se degrada: explota.
ne = 16;
[nodeE,elemE] = squaremesh_basico(ne);
dtcE = dtCritico(nodeE,elemE,0);
he   = 1/ne;
ueE  = uex(nodeE(:,1),nodeE(:,2),T);
fprintf('\n===== ESTABILIDAD  (n = %d, h = %.4f) =====\n', ne, he);
fprintf('dt* = 2/lambda_max = %.3e   (h^2/12.93 = %.3e, condensada h^2/4 = %.3e)\n', ...
    dtcE, he^2/12.93, he^2/4);
fprintf('   pasos      dt      dt/dt*    err atras     err CN     err adelante\n');
for ns = [50 100 200 1000]
    dtk = T/ns;
    e1 = norm(resolverBE(nodeE,elemE,dtk,T,u0f,f) - ueE, inf);
    e2 = norm(resolverCN(nodeE,elemE,dtk,T,u0f,f) - ueE, inf);
    e3 = norm(resolverFE(nodeE,elemE,dtk,T,u0f,f) - ueE, inf);
    fprintf('%8d  %.3e  %6.2f   %.3e   %.3e   %.3e\n', ns, dtk, dtk/dtcE, e1, e2, e3);
end

%% --- 7. Paso enorme ---------------------------------------------------
%  Un solo paso que cubre todo el intervalo. El explicito explotaria sin
%  remedio; los otros dos siguen acotados, apenas imprecisos.
fprintf('\n===== PASO ENORME  (n = %d) =====\n', n);
fprintf('   pasos      dt      err atras    err CN\n');
for ns = [1 2 5 200]
    e1 = norm(resolverBE(node,elem,T/ns,T,u0f,f) - ue, inf);
    e2 = norm(resolverCN(node,elem,T/ns,T,u0f,f) - ue, inf);
    fprintf('%8d  %.3e  %.3e  %.3e\n', ns, T/ns, e1, e2);
end
fprintf('En esta malla el explicito pediria dt <= %.2e, o sea al menos\n', dtc);
fprintf('%d pasos para llegar a T = %g.\n', ceil(T/dtc), T);

%% --- 8. Fuente dependiente del tiempo ---------------------------------
%  Comprobacion de que el promedio del lado derecho conserva el orden 2. Se
%  fabrica una solucion que no es la de decaimiento natural,
%
%      u = e^{-t} sin(pi x) sin(pi y)   =>   f = (2 pi^2 - 1) e^{-t} sin sin
%
%  de modo que f si dependa del tiempo. Si el lado derecho se evaluara solo en
%  t^{n+1}, en vez de promediarlo, Crank-Nicolson bajaria a orden 1.
fprintf('\n===== FUENTE DEPENDIENTE DEL TIEMPO  (n = %d) =====\n', n);
ft   = @(x,y,t) (2*pi^2-1)*exp(-t).*sin(pi*x).*sin(pi*y);
u0t  = @(x,y) sin(pi*x).*sin(pi*y);
reft = resolverCN(node,elem,T/40000,T,u0t,ft);     % referencia, misma malla
fprintf('   pasos      dt       err CN     orden\n');
ep = NaN;
for ns = [10 20 40 80]
    e = norm(resolverCN(node,elem,T/ns,T,u0t,ft) - reft, inf);
    fprintf('%8d  %.2e  %.3e  %7.2f\n', ns, T/ns, e, log2(ep/e));
    ep = e;
end
uexT = @(x,y) exp(-T).*sin(pi*x).*sin(pi*y);
fprintf('error contra la solucion exacta con dt = T/200: %.4e\n', ...
    norm(resolverCN(node,elem,T/200,T,u0t,ft) - uexT(node(:,1),node(:,2)), inf));

%% ====================== funciones locales ============================

function [node,elem] = squaremesh_basico(M)
%SQUAREMESH_BASICO Malla estructurada de (0,1)^2, M subintervalos por lado.
%   node: (M+1)^2 x 2   elem: 2*M^2 x 3, orientados en sentido antihorario.
[X,Y] = meshgrid(linspace(0,1,M+1));
node  = [X(:), Y(:)];   % el nodo (x_i,y_j) ocupa la fila (i-1)*(M+1)+j
elem  = zeros(2*M^2,3);
k = 0;
for i = 1:M
    for j = 1:M
        n1 = (i-1)*(M+1) + j;
        n2 = n1 + (M+1);
        n3 = n2 + 1;
        n4 = n1 + 1;
        elem(k+1,:) = [n1 n2 n3];
        elem(k+2,:) = [n1 n3 n4];
        k = k + 2;
    end
end
end

function free = nodosLibres(elem,NV)
%NODOSLIBRES  Nodos que no estan sobre la frontera.
%   Una arista es de frontera si pertenece a un solo triangulo.
totalEdge  = sort([elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])], 2);
[edge,~,j] = unique(totalEdge,'rows');
conteo     = accumarray(j,1);
bdNode     = unique(edge(conteo==1,:));
free       = setdiff(1:NV, bdNode);
end

function [A,M,b,area] = ensamblar(node,elem,f)
%ENSAMBLAR  Rigidez, masa y lado derecho P1 en 2D. Es assembling.m.
N = size(node,1); NT = size(elem,1);
ii = zeros(9*NT,1); jj = zeros(9*NT,1); sA = zeros(9*NT,1); sM = zeros(9*NT,1);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);   % lado opuesto al vertice 3
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
index = 0;
for i = 1:3
    for j = 1:3
        ii(index+1:index+NT) = elem(:,i);
        jj(index+1:index+NT) = elem(:,j);
        sA(index+1:index+NT) = dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        sM(index+1:index+NT) = area*((i==j)+1)/12;
        index = index + NT;
    end
end
A = sparse(ii,jj,sA,N,N);
M = sparse(ii,jj,sM,N,N);
barycenter = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
bt = f(barycenter(:,1),barycenter(:,2))/3;
bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[N 1]);
end

function u = resolverTheta(node,elem,dt,T,u0f,f,theta)
%RESOLVERTHETA  Metodo theta:
%
%   (M + theta dt A) u^{n+1}
%       = (M - (1-theta) dt A) u^n + dt [theta b^{n+1} + (1-theta) b^n]
%
%   theta = 0 explicito, 1 implicito, 1/2 Crank-Nicolson.
%
%   La fuente puede darse de dos formas y se detecta sola:
%       f = @(x,y)   -> estacionaria
%       f = @(x,y,t) -> dependiente del tiempo
%
%   EL LADO DERECHO LLEVA LA MISMA COMBINACION QUE A*u, y no es un capricho.
%   Integrando el sistema sobre [t_n, t_{n+1}] queda
%
%       M (u^{n+1}-u^n) + A int u dt = int b dt
%
%   y el metodo theta es elegir una cuadratura para esas dos integrales. Con
%   theta = 1/2 la de la izquierda es la regla del trapecio, de error local
%   O(dt^3). Si la de la derecha se aproximara por el rectangulo dt*b^{n+1},
%   que es solo O(dt^2), mandaria el termino mas debil y Crank-Nicolson caeria
%   a orden 1. Medido con una fuente dependiente del tiempo, el error llega a
%   ser cuatro ordenes de magnitud mayor.
%
%   Con f estacionaria b^{n+1} = b^n y la combinacion se reduce a b, que es el
%   caso de este archivo.
NV = size(node,1);
free = nodosLibres(elem,NV);

% nargin de un handle anonimo dice cuantos argumentos declara
variaEnTiempo = (nargin(f) == 3);
if variaEnTiempo
    [A,M] = ensamblar(node,elem,@(x,y) zeros(size(x)));
else
    [A,M,b] = ensamblar(node,elem,f);
end

u = u0f(node(:,1),node(:,2));
u(setdiff(1:NV,free)) = 0;            % Dirichlet homogeneo
if T == 0, return, end

nsteps = max(1, round(T/dt));
dt = T/nsteps;                        % ajusta dt para caer exactamente en T

% Ninguna de las dos matrices depende del paso: se arman una sola vez.
Mf = M(free,free);
Af = A(free,free);
dK = decomposition(Mf + theta*dt*Af, 'chol');   % simetrica def. positiva
R  = Mf - (1-theta)*dt*Af;

if variaEnTiempo
    bNuevo = ladoDerecho(node,elem,@(x,y) f(x,y,0));
end

for k = 1:nsteps
    if variaEnTiempo
        % el b^{n+1} del paso anterior es el b^n de este: se reaprovecha
        bViejo = bNuevo;
        bNuevo = ladoDerecho(node,elem,@(x,y) f(x,y,k*dt));
        bk = theta*bNuevo + (1-theta)*bViejo;
    else
        bk = b;
    end
    u(free) = dK \ (R*u(free) + dt*bk(free));
end
end

function b = ladoDerecho(node,elem,f)
%LADODERECHO  Solo el vector b del ensamblaje, para rehacerlo en cada paso
%   cuando la fuente depende del tiempo. Es la parte final de ensamblar.
N = size(node,1);
v1 = node(elem(:,2),:)-node(elem(:,1),:);
v2 = node(elem(:,3),:)-node(elem(:,1),:);
area = 0.5*abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1));
barycenter = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
bt = f(barycenter(:,1),barycenter(:,2))/3;
bt = bt.*repmat(area,1,3);
b = accumarray(elem(:),bt(:),[N 1]);
end

function u = resolverBE(node,elem,dt,T,u0f,f)
%RESOLVERBE  Euler hacia atras (theta = 1). Incondicionalmente estable.
u = resolverTheta(node,elem,dt,T,u0f,f,1.0);
end

function u = resolverCN(node,elem,dt,T,u0f,f)
%RESOLVERCN  Crank-Nicolson (theta = 1/2). Orden 2 en el tiempo.
u = resolverTheta(node,elem,dt,T,u0f,f,0.5);
end

function u = resolverFE(node,elem,dt,T,u0f,f)
%RESOLVERFE  Euler hacia adelante (theta = 0). Estable solo si dt <= dtCritico.
u = resolverTheta(node,elem,dt,T,u0f,f,0.0);
end

function lh = lambdaMin(node,elem)
%LAMBDAMIN  Menor valor propio de A v = lambda M v.
%   Es la tasa de decaimiento del modo mas lento de la solucion discreta, el
%   analogo de 2*pi^2 en el continuo. Tiende a 2*pi^2 + O(h^2) por arriba.
NV = size(node,1);
[A,M] = ensamblar(node,elem,@(x,y) zeros(size(x)));
free = nodosLibres(elem,NV);
lh = eigs(A(free,free), M(free,free), 1, 'smallestabs');
end

function dtc = dtCritico(node,elem,theta)
%DTCRITICO  Paso maximo admisible del metodo theta.
%   De |g| <= 1 sale dt*lambda_max*(1-2 theta) <= 2. Para theta >= 1/2 no hay
%   restriccion. lambda_max es el mayor valor propio de A v = lambda M v; en
%   esta malla tiende a 25.856/h^2, de donde dtc -> h^2/12.93 cuando theta = 0.
if theta >= 0.5
    dtc = Inf;
    return
end
NV = size(node,1);
[A,M] = ensamblar(node,elem,@(x,y) zeros(size(x)));
free = nodosLibres(elem,NV);
lmax = eigs(A(free,free), M(free,free), 1, 'largestabs');
dtc = 2/((1-2*theta)*lmax);
end

function e = errorL2(node,elem,u,uexT)
%ERRORL2  Norma L2 del error, con la regla de grado 5 (7 puntos) del capitulo.
a1=0.470142064105115; b1=0.059715871789770;
a2=0.101286507323456; b2=0.797426985353087;
L = [1/3 1/3 1/3; b1 a1 a1; a1 b1 a1; a1 a1 b1; b2 a2 a2; a2 b2 a2; a2 a2 b2];
W = [9/40; repmat(0.132394152788506,3,1); repmat(0.125939180544827,3,1)];

[~,~,~,area] = ensamblar(node,elem,@(x,y) zeros(size(x)));
e2 = 0;
for q = 1:size(L,1)
    % punto de cuadratura y valor de u_h ahi, en coordenadas baricentricas
    xq = L(q,1)*node(elem(:,1),1) + L(q,2)*node(elem(:,2),1) + L(q,3)*node(elem(:,3),1);
    yq = L(q,1)*node(elem(:,1),2) + L(q,2)*node(elem(:,2),2) + L(q,3)*node(elem(:,3),2);
    uq = L(q,1)*u(elem(:,1))      + L(q,2)*u(elem(:,2))      + L(q,3)*u(elem(:,3));
    e2 = e2 + sum( W(q) * (uq - uexT(xq,yq)).^2 .* area );
end
e = sqrt(e2);
end
