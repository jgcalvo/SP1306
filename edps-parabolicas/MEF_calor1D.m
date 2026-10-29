%MEF_CALOR1D  Ecuacion del calor en 1D con elementos P1: Euler implicito y explicito.
%
%   Resuelve
%
%       u_t = u_xx            en (0,1) x (0,T]
%       u(0,t) = u(1,t) = 0
%       u(x,0) = u0(x)
%
%   con elementos finitos P1 en el espacio, avanzando en el tiempo de dos
%   maneras: Euler hacia atras (implicito) y Euler hacia adelante (explicito).
%
%   FORMULACION DEBIL.  Multiplicando por v que se anula en los extremos e
%   integrando por partes solo el termino de segundo orden,
%
%       (u_t, v) + (u_x, v_x) = (f, v)      para toda v
%
%   Al escribir u_h(x,t) = sum_j u_j(t) phi_j(x) aparecen las DOS matrices del
%   ensamblaje: la de masa viene del termino temporal y la de rigidez del
%   termino espacial,
%
%       M u'(t) + A u(t) = b,     M_ij = int phi_i phi_j,  A_ij = int phi_i' phi_j'
%
%   Esta es la diferencia con el problema estacionario, donde M no se usa.
%
%   LOS DOS ESQUEMAS.  Segun donde se evalue el termino A u:
%
%       hacia atras:    M (u^{n+1}-u^n)/dt + A u^{n+1} = b
%                       (M + dt A) u^{n+1} = M u^n + dt b
%
%       hacia adelante: M (u^{n+1}-u^n)/dt + A u^n     = b
%                       M u^{n+1} = (M - dt A) u^n + dt b
%
%   Un detalle que conviene no pasar por alto: en elementos finitos el metodo
%   "explicito" TAMBIEN pide resolver un sistema en cada paso, porque la matriz
%   de masa M no es diagonal. Lo explicito es el tratamiento del operador
%   espacial, no el costo del paso.
%
%   En los dos casos la matriz que hay que invertir -- M + dt A en uno, M en el
%   otro -- NO cambia de un paso al otro, porque ni la malla ni dt dependen del
%   tiempo. Se factoriza UNA sola vez, fuera del ciclo, y cada paso es solo una
%   sustitucion. Es el mismo criterio que lleva a ensamblar de una sola vez en
%   vez de entrada por entrada.
%
%   ESTABILIDAD.  Aqui esta la diferencia de fondo. Si lambda_h son los valores
%   propios de A v = lambda M v, el factor de amplificacion de cada modo es
%
%       hacia atras:     1/(1 + dt lambda)   < 1 siempre
%       hacia adelante:  1 - dt lambda       en valor absoluto <= 1 solo si
%                                            dt <= 2/lambda_max
%
%   El implicito es incondicionalmente estable: dt y h se eligen por precision.
%   El explicito no. Para P1 en una malla uniforme se puede calcular
%
%       lambda_max  ->  12/h^2      de donde      dt* = 2/lambda_max = h^2/6
%
%   ENSAMBLAJE.  La funcion local ensamblar1D es la version en una dimension de
%   assembling.m: misma estructura de tripletes y mismo ciclo sobre los pares de
%   funciones de forma, con 2x2 = 4 entradas por elemento en vez de 3x3 = 9.
%
%   Sobre un elemento de longitud h las matrices locales son
%
%       A_loc = (1/h) [ 1 -1 ; -1  1 ]       M_loc = (h/6) [ 2 1 ; 1 2 ]
%
%   que en forma compacta, para el ciclo vectorizado, se escriben
%
%       sA = (2*(i==j)-1)./h            sM = h.*((i==j)+1)/6
%
%   La segunda es la hermana en 1D de la formula area*((i==j)+1)/12 del codigo
%   en 2D.
%
%   ORDENES.  Los dos esquemas son de orden 1 en dt y orden 2 en h. Para este
%   dato inicial el error espacial se puede rastrear hasta su origen: el valor
%   propio discreto mas pequeno cumple
%
%       lambda_h = pi^2 + pi^4 h^2/12 + O(h^4)
%
%   asi que la solucion numerica decae un poco mas rapido de lo debido, y el
%   error en el instante t vale aproximadamente t*(lambda_h - pi^2)*e^{-pi^2 t}.
%   Con h = 0.1 y t = 0.1 eso da 3.0e-3, que es justo lo que mide la seccion 4.

clear; clc

%% --- 1. Datos del problema -------------------------------------------
T    = 0.1;                            % tiempo final
u0f  = @(x) sin(pi*x);                 % dato inicial
f    = @(x) zeros(size(x));            % sin fuente: u_t = u_xx
uex  = @(x,t) exp(-pi^2*t).*sin(pi*x); % solucion exacta correspondiente

% Con u0 = sin(pi x) la solucion exacta es e^{-pi^2 t} sin(pi x): la forma se
% conserva y solo decae la amplitud.

N  = 101;                              % nodos
dt = T/100;                            % paso de tiempo

[node,elem] = malla1D(0,1,N);
uBE = resolverBE(node,elem,dt,T,u0f,f);

fprintf('N = %d nodos, dt = %.2e, %d pasos\n', N, dt, round(T/dt));
fprintf('dt critico del explicito en esta malla: %.3e\n', dtCritico(node,elem));
fprintf('\nEuler hacia atras:\n');
fprintf('  error max en los nodos, t = T: %.4e\n', norm(uBE - uex(node,T), inf));
fprintf('  error en L2,            t = T: %.4e\n', ...
    errorL2(node,elem,uBE,@(x) uex(x,T)));

%% --- 2. Evolucion en el tiempo ----------------------------------------
figure
instantes = [0 0.005 0.02 0.05 0.1];
colores = lines(numel(instantes));
hold on
for k = 1:numel(instantes)
    tk = instantes(k);
    uk = resolverBE(node,elem,dt,tk,u0f,f);
    plot(node, uk, '-',  'Color', colores(k,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('t = %g', tk));
    plot(node, uex(node,tk), '--k', 'HandleVisibility','off');
end
hold off
xlabel('x'); ylabel('u'); legend show
title('Euler implicito (lineas) vs solucion exacta (punteado)')

%% --- 3. Convergencia en el tiempo -------------------------------------
%  Se fija una malla fina para que el error espacial no contamine la medicion.
%  Solo el implicito: en esta malla el explicito necesitaria dt < 4e-8.
fprintf('\n===== ORDEN EN EL TIEMPO, IMPLICITO  (N = 2001 fijo) =====\n');
fprintf('   pasos      dt        err_max     orden\n');
[nodeF,elemF] = malla1D(0,1,2001);
ep = NaN;
for ns = [10 20 40 80 160]
    dtk = T/ns;
    uk  = resolverBE(nodeF,elemF,dtk,T,u0f,f);
    e   = norm(uk - uex(nodeF,T), inf);
    fprintf('%8d  %.3e  %.4e  %7.2f\n', ns, dtk, e, log2(ep/e));
    ep = e;
end
fprintf('Euler hacia atras es de orden 1: el error se parte a la mitad.\n');

%% --- 4. Convergencia en el espacio ------------------------------------
%  Ahora al reves: dt muy pequeno para aislar el error de la discretizacion
%  espacial. Hacen falta muchisimos pasos justamente porque Euler hacia atras
%  es de orden 1: para que su error quede por debajo del espacial hay que
%  exagerar. Con dt mas grande los ordenes de abajo se inflan por encima de 2,
%  que es la senal de que las dos fuentes de error ya son comparables.
fprintf('\n===== ORDEN EN EL ESPACIO, IMPLICITO  (dt = T/200000 fijo) =====\n');
fprintf('    N        h        err_max    orden      err_L2    orden\n');
dtf = T/200000;
emp = NaN; e2p = NaN;
for Nk = [11 21 41 81]
    [nk,ek] = malla1D(0,1,Nk);
    uk = resolverBE(nk,ek,dtf,T,u0f,f);
    em = norm(uk - uex(nk,T), inf);
    e2 = errorL2(nk,ek,uk,@(x) uex(x,T));
    fprintf('%5d  %.3e  %.3e  %7.2f   %.3e  %7.2f\n', ...
        Nk, 1/(Nk-1), em, log2(emp/em), e2, log2(e2p/e2));
    emp = em; e2p = e2;
end

%% --- 5. Implicito contra explicito: la restriccion de estabilidad -----
%  Misma malla y mismo dt para los dos. El implicito funciona en toda la
%  columna; el explicito solo cuando dt queda por debajo de dt*, y apenas se
%  pasa no se degrada: explota.
Ne = 41;
[nodeE,elemE] = malla1D(0,1,Ne);
dtc = dtCritico(nodeE,elemE);
he  = 1/(Ne-1);
fprintf('\n===== IMPLICITO vs EXPLICITO  (N = %d, h = %.4f) =====\n', Ne, he);
fprintf('dt* = 2/lambda_max = %.3e   (h^2/6 = %.3e, h^2/2 = %.3e)\n', ...
    dtc, he^2/6, he^2/2);
fprintf('   pasos      dt      dt/dt*     err implic.   err explic.\n');
for ns = [100 500 1000 2000]
    dtk = T/ns;
    eI = norm(resolverBE(nodeE,elemE,dtk,T,u0f,f) - uex(nodeE,T), inf);
    eE = norm(resolverFE(nodeE,elemE,dtk,T,u0f,f) - uex(nodeE,T), inf);
    fprintf('%8d  %.3e  %7.2f   %10.3e    %10.3e\n', ns, dtk, dtk/dtc, eI, eE);
end
% El error del implicito no baja monotonamente en esa columna, y no es un
% error del codigo: Euler hacia atras decae de menos (1/(1+dt*lam) > e^{-dt*lam})
% mientras que la discretizacion espacial decae de mas (lambda_h > pi^2). Los
% dos errores tienen signos opuestos y se cancelan cerca de dt = 1e-4.
fprintf('El error implicito no es monotono: el error temporal y el espacial\n');
fprintf('tienen signos opuestos y se cancelan cerca de dt = 1e-4.\n');
fprintf('Donde el explicito es estable, los dos dan practicamente lo mismo:\n');
fprintf('son el mismo metodo de orden 1, y no hay razon para pagar la\n');
fprintf('restriccion sobre dt.\n');

%% --- 6. Estabilidad incondicional del implicito -----------------------
%  Un solo paso que cubre todo el intervalo de tiempo. El explicito explotaria
%  sin remedio; aca la solucion sigue siendo razonable, apenas imprecisa, que
%  es exactamente lo que significa "estable pero inexacto".
fprintf('\n===== IMPLICITO CON dt ENORME  (N = %d) =====\n', N);
fprintf('   pasos      dt       err_max    max|u|\n');
for ns = [1 2 5 100]
    uk = resolverBE(node,elem,T/ns,T,u0f,f);
    fprintf('%8d  %.3e  %.3e  %.4f\n', ns, T/ns, ...
        norm(uk-uex(node,T),inf), max(abs(uk)));
end
lim = dtCritico(node,elem);
fprintf('En esta malla el explicito pediria dt <= %.2e, o sea al menos\n', lim);
fprintf('%d pasos para llegar a T = %g.\n', ceil(T/lim), T);

%% ====================== funciones locales ============================

function [node,elem] = malla1D(a,b,N)
%MALLA1D  Malla de N nodos en [a,b]: node es Nx1 y elem es (N-1)x2.
node = linspace(a,b,N)';
elem = [(1:N-1)', (2:N)'];
end

function [A,M,b,h] = ensamblar1D(node,elem,f)
%ENSAMBLAR1D  Rigidez, masa y lado derecho P1 en 1D.
%   Version en una dimension de assembling.m: mismo formato de tripletes y
%   mismo ciclo sobre los pares (i,j), con 4 entradas por elemento.
N = size(node,1); NE = size(elem,1);
ii = zeros(4*NE,1); jj = zeros(4*NE,1); sA = zeros(4*NE,1); sM = zeros(4*NE,1);
h = node(elem(:,2)) - node(elem(:,1));      % longitud de cada elemento
index = 0;
for i = 1:2
    for j = 1:2
        ii(index+1:index+NE) = elem(:,i);
        jj(index+1:index+NE) = elem(:,j);
        sA(index+1:index+NE) = (2*(i==j)-1)./h;   % (1/h)[1 -1; -1 1]
        sM(index+1:index+NE) = h*((i==j)+1)/6;    % (h/6)[2 1; 1 2]
        index = index + NE;
    end
end
A = sparse(ii,jj,sA,N,N);
M = sparse(ii,jj,sM,N,N);
% lado derecho con la regla del punto medio, exacta para grado 1
medio = (node(elem(:,1)) + node(elem(:,2)))/2;
bt = f(medio)/2;
bt = bt.*repmat(h,1,2);
b = accumarray(elem(:),bt(:),[N 1]);
end

function u = resolverBE(node,elem,dt,T,u0f,f)
%RESOLVERBE  Euler hacia atras:  (M + dt A) u^{n+1} = M u^n + dt b.
N = size(node,1);
[A,M,b,~] = ensamblar1D(node,elem,f);

free = 2:N-1;                         % Dirichlet homogeneo en ambos extremos
u = u0f(node);
u([1 N]) = 0;
if T == 0, return, end

nsteps = max(1, round(T/dt));
dt = T/nsteps;                        % ajusta dt para caer exactamente en T

% La matriz del sistema no depende del paso: se factoriza una sola vez.
Mf = M(free,free);
dK = decomposition(Mf + dt*A(free,free), 'chol');   % simetrica def. positiva

for n = 1:nsteps
    u(free) = dK \ (Mf*u(free) + dt*b(free));
end
end

function u = resolverFE(node,elem,dt,T,u0f,f)
%RESOLVERFE  Euler hacia adelante:  M u^{n+1} = (M - dt A) u^n + dt b.
%   Aunque el esquema sea explicito, con masa consistente hay que resolver con
%   M en cada paso. M tampoco cambia, asi que se factoriza una sola vez.
%   Es estable solo si dt <= dtCritico(node,elem).
N = size(node,1);
[A,M,b,~] = ensamblar1D(node,elem,f);

free = 2:N-1;
u = u0f(node);
u([1 N]) = 0;
if T == 0, return, end

nsteps = max(1, round(T/dt));
dt = T/nsteps;

Mf = M(free,free);
Af = A(free,free);
dM = decomposition(Mf, 'chol');

for n = 1:nsteps
    u(free) = dM \ (Mf*u(free) - dt*(Af*u(free)) + dt*b(free));
end
end

function dtc = dtCritico(node,elem)
%DTCRITICO  Paso maximo que admite Euler hacia adelante: 2/lambda_max.
%   lambda_max es el mayor valor propio de A v = lambda M v. En una malla
%   uniforme tiende a 12/h^2, de donde dtc -> h^2/6.
N = size(node,1);
[A,M] = ensamblar1D(node,elem,@(x) zeros(size(x)));
free = 2:N-1;
lmax = eigs(A(free,free), M(free,free), 1, 'largestabs');
dtc = 2/lmax;
end

function e = errorL2(node,elem,u,uexT)
%ERRORL2  Norma L2 del error, con Gauss de 3 puntos en cada elemento.
%   Exacta para grado 5, de sobra para comparar una funcion lineal a trozos
%   contra una suave.
xi = [-sqrt(3/5); 0; sqrt(3/5)];
w  = [5/9; 8/9; 5/9];

x1 = node(elem(:,1)); x2 = node(elem(:,2));
u1 = u(elem(:,1));    u2 = u(elem(:,2));
h  = x2 - x1;

e2 = 0;
for q = 1:numel(xi)
    xq = (x1+x2)/2 + h/2*xi(q);              % punto de cuadratura
    uq = u1.*(1-xi(q))/2 + u2.*(1+xi(q))/2;  % u_h es lineal en el elemento
    e2 = e2 + sum( w(q) * (uq - uexT(xq)).^2 .* h/2 );
end
e = sqrt(e2);
end
