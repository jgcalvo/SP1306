% MDF v3 - Estudio de convergencia, sin ciclos
%
% Mismo problema y mismo estudio que MDF_v1_converg.m y MDF_v2_converg.m, para
% el caso de condiciones de frontera homogeneas:
%
%     -Lap(u) = f    en  Omega = (0,1)x(0,1)
%           u = 0    sobre la frontera de Omega
%
% Las versiones anteriores recorren los nodos uno por uno y van llenando las
% listas ii, jj, ss. Aca no hay ciclo sobre los nodos: ni para la matriz, ni
% para el lado derecho.
%
%
% LA MATRIZ COMO PRODUCTO DE KRONECKER
%
% En una dimension, la formula de segundas diferencias sobre los n = M-1 nodos
% interiores da la matriz tridiagonal
%
%     T = tridiag(-1, 2, -1)   de tamano n x n
%
% El laplaciano en 2D es la suma de la segunda derivada en x y la segunda
% derivada en y, y esa separacion se hereda a la matriz. La numeracion
%
%     pos = (j-1)*(M-1) + i
%
% recorre primero i (la direccion x) y despues j, asi que el vector de
% incognitas queda partido en n bloques, uno por cada valor de j:
%
%     u = [ u(:,1) ; u(:,2) ; ... ; u(:,n) ]
%
% Con esa estructura por bloques:
%
%   - la derivada en x actua DENTRO de cada bloque, sin mezclarlos: es una
%     matriz diagonal por bloques con T en cada bloque, o sea kron(I,T)
%
%   - la derivada en y conecta el bloque j con los bloques j-1 y j+1, nodo con
%     nodo (mismo i): los bloques diagonales son 2*I y los vecinos -I, o sea
%     kron(T,I)
%
% Sumando,
%
%     A = kron(I,T) + kron(T,I)
%
% La diagonal queda 2+2 = 4 y cada vecino -1, que es exactamente la formula de
% 5 puntos. Es la misma matriz que arman las versiones v1 y v2, entrada por
% entrada, pero construida de una sola vez.
%
% Importante: I y T se construyen RALAS (speye, spdiags). kron de dos matrices
% ralas devuelve una rala; si se usara eye(n) densa, kron intentaria reservar
% una matriz de (M-1)^2 x (M-1)^2 elementos y no cabria en memoria.
%
%
% EL LADO DERECHO SIN CICLO
%
% Hay que llenar b(pos) = h^2 * f(x_i, y_j) respetando el mismo orden. Como
% b = B(:) apila la matriz B por columnas, y el indice de columna es j, basta
% construir B con i en las filas y j en las columnas:
%
%     B(i,j) = h^2 * f(x_i, y_j)
%
% Eso es justo lo que produce ndgrid, no meshgrid: [XX,YY] = ndgrid(x,y) da
% XX(i,j) = x(i) e YY(i,j) = y(j). meshgrid devuelve la transpuesta de eso y
% dejaria el vector en el orden equivocado.

% entradas
uex = @(x,y) sin(pi*x).*sin(pi*y);        % sol. exacta
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y); % lado derecho

err = [];                                 % error de cada malla
hh = [];                                  % paso h de cada malla

for M = 2.^(2:7)                          % cantidad intervalos
    tic
    h = 1/M;                              % paso de la malla
    n = M-1;                              % nodos interiores por direccion

    % --- matriz: laplaciano 1D y productos de Kronecker
    e = ones(n,1);
    T = spdiags([-e 2*e -e], -1:1, n, n); % tridiag(-1,2,-1), rala
    I = speye(n);
    A = kron(I,T) + kron(T,I);            % formula de 5 puntos, de una sola vez

    % --- lado derecho y solucion exacta, sin ciclo
    x = (1:n)*h;                          % coordenadas de los nodos interiores
    [XX,YY] = ndgrid(x,x);                % XX(i,j) = x_i,  YY(i,j) = y_j
    B      = h^2*f(XX,YY);                % B(i,j) = h^2 f(x_i,y_j)
    UexMat = uex(XX,YY);                  % referencia para medir el error
    b = B(:);                             % apila por columnas: orden pos

    uh = A\b;                             % resuelve el sistema
    toc

    hh = [hh, h];
    err = [err, norm(uh-UexMat(:),inf)];
end

%%
% en log-log el error debe quedar paralelo a h^2, es decir con pendiente 2
loglog(hh,err,'r'), hold on
loglog(hh,hh.^2,'k')
legend('err','h^2')
