% MDF - Metodo de diferencias finitas para la ecuacion de Poisson en 2D
%
% Resuelve el problema con condiciones de frontera homogeneas
%
%     -Lap(u) = f    en  Omega = (0,1)x(0,1)
%           u = 0    sobre la frontera de Omega
%
% Discretizacion: malla uniforme de M intervalos por direccion, h = 1/M, con
% nodos (x_i,y_j) = (i*h, j*h). La formula de 5 puntos aproxima el laplaciano
% en el nodo (i,j) como
%
%     -Lap(u) ~ ( 4*u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1) ) / h^2
%
% Multiplicando la ecuacion por h^2 se evita dividir, y el sistema queda
%
%     4*u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1) = h^2 * f(i,j)
%
% Las incognitas son solo los nodos interiores, i,j = 1..M-1: los de la
% frontera valen 0 y por eso no entran al sistema. En total dimA = (M-1)^2.
%
% Numeracion de los nodos: a cada par (i,j) le corresponde una unica fila
%
%     pos = (j-1)*(M-1) + i
%
% Con esta numeracion los vecinos en x quedan en pos-1 y pos+1 (contiguos), 
% y los vecinos en y quedan a distancia M-1. 
%
% Ensamblaje: en vez de escribir A(fila,columna) = valor una entrada a la 
% vez, se van acumulando tres listas paralelas ii (filas), jj (columnas) y 
% ss (valores). Al final sparse(ii,jj,ss,...) construye la matriz de una 
% sola vez. Cada terna (ii(k), jj(k), ss(k)) dice
% "en la fila ii(k), columna jj(k), va el valor ss(k)".

% entradas
uex = @(x,y) sin(2*pi*x).*sin(pi*y);        % sol. exacta
f   = @(x,y) 5*pi^2*sin(2*pi*x).*sin(pi*y); % lado derecho


M = 32;                                     % cantidad intervalos
h = 1/M;                                    % paso de la malla
dimA = (M-1)^2;                             % cantidad de nodos interiores
% vectores para generar A
ii = []; jj = []; ss = [];
b = zeros(dimA,1);                          % lado derecho del sistema
UexMat = zeros(M-1,M-1);                        % sol. exacta evaluada en la malla

% recorrido de los nodos interiores: j avanza en y, i avanza en x
for j = 1:M-1
    for i = 1:M-1
        pos = (j-1)*(M-1)+i;                % fila del sistema para el nodo (i,j)
        b(pos) = h^2*f(i*h,j*h);            % lado derecho, ya multiplicado por h^2
        UexMat(i,j) = uex(i*h,j*h);             % referencia para medir el error

        % termino diagonal: el propio nodo (i,j)
        ii = [ii, pos];
        jj = [jj, pos];
        ss = [ss, 4];

        % vecino izquierdo (i-1,j)
        if(i~=1)
            ii = [ii, pos-1];
            jj = [jj, pos];
            ss = [ss, -1];
        else
            % si i==1 el vecino esta sobre la frontera x=0, donde u=0: su
            % termino se anula y no aporta nada al sistema. Con condiciones
            % NO homogeneas se debe restar el valor de frontera a b(pos)
        end

        % vecino derecho (i+1,j); si i==M-1 cae en la frontera x=1
        if(i~=M-1)
            ii = [ii, pos+1];
            jj = [jj, pos];
            ss = [ss, -1];
        end

        % vecino inferior (i,j-1); si j==1 cae en la frontera y=0
        if(j~=1)
            ii = [ii, pos];
            jj = [jj, pos-(M-1)];
            ss = [ss, -1];
        end

        % vecino superior (i,j+1); si j==M-1 cae en la frontera y=1
        if(j~=M-1)
            ii = [ii, pos];
            jj = [jj, pos+(M-1)];
            ss = [ss, -1];
        end
    end
end

A = sparse(ii,jj,ss,dimA,dimA);             % arma la matriz de una sola vez
%Afull = full(A);
uh = A\b;                                   % resuelve el sistema

% visualización
xx = linspace(0,1,M+1);                     % malla completa, con frontera
[X,Y] = meshgrid(xx,xx);
Z = zeros(M+1,M+1);                         % los ceros del borde son la cond. de frontera
% reshape devuelve la solucion a forma de matriz; la transpuesta es necesaria
% porque meshgrid indexa Z(fila,columna) = Z(y,x), al reves que UexMat(i,j) = u(x,y)
Z(2:end-1,2:end-1) = reshape(uh,M-1,M-1)';
surf(X,Y,Z)
xlabel('x')
ylabel('y')

% error norma inf
% restar uex - uh
% (ojo: tal como esta, norm(v) es la norma 2 del vector; para la norma
%  infinito hace falta norm(v,inf), como en MDF_v1_converg.m)
norm(uh-UexMat(:))
