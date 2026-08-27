% MDF v2 - Diferencias finitas en 2D, con memoria preasignada
%
% Misma discretizacion que MDF_v1_test.m: formula de 5 puntos para
%
%     -Lap(u) = f    en  Omega = (0,1)x(0,1)
%           u = 0    sobre la frontera de Omega
%
% y la misma numeracion de nodos, pos = (j-1)*(M-1) + i, con i el indice rapido.
% La matriz que se obtiene es exactamente la misma que en la v1.
%
% QUE CAMBIA RESPECTO A LA v1
%
% En la v1 los vectores crecen de a un elemento: ii = [ii, pos]. Cada una de
% esas asignaciones obliga a MATLAB a pedir un bloque de memoria mas grande y a
% copiar todo lo que ya habia. Si al final hay N entradas, el costo total de
% armar los vectores es O(N^2) en vez de O(N), y eso domina el tiempo del
% programa: para M=128 tarda mas en armar las listas que en resolver el sistema.
%
% Aca se reserva de entrada el espacio exacto y se llena con un contador k. La
% cuenta de cuanto espacio hace falta:
%
%   - cada uno de los (M-1)^2 nodos interiores aporta 5 entradas (el propio
%     nodo y sus 4 vecinos)                                 -> 5*(M-1)^2
%   - pero los nodos pegados a la frontera pierden un vecino cada uno. Hay M-1
%     nodos con i=1, otros M-1 con i=M-1, otros M-1 con j=1 y otros M-1 con
%     j=M-1                                                 -> -4*(M-1)
%
%     nnz = 5*(M-1)^2 - 4*(M-1)
%
% (para M=32 esto da 4681, que es exactamente nnz(A))
%
% Este conteo vale para condiciones de frontera HOMOGENEAS, que es el caso de
% este archivo: los vecinos que caen sobre la frontera valen 0 y simplemente no
% generan entrada.

% entradas
uex = @(x,y) sin(2*pi*x).*sin(pi*y);        % sol. exacta
f   = @(x,y) 5*pi^2*sin(2*pi*x).*sin(pi*y); % lado derecho

M = 32;                                     % cantidad intervalos
h = 1/M;                                    % paso de la malla
dimA = (M-1)^2;                             % cantidad de nodos interiores

% espacio exacto para las tres listas del formato triplete
nnzA = 5*(M-1)^2 - 4*(M-1);
ii = zeros(nnzA,1);                         % filas
jj = zeros(nnzA,1);                         % columnas
ss = zeros(nnzA,1);                         % valores
k  = 0;                                     % cuantas entradas se han escrito

b  = zeros(dimA,1);                         % lado derecho del sistema
UexMat = zeros(M-1,M-1);                        % sol. exacta evaluada en la malla

% recorrido de los nodos interiores: j avanza en y, i avanza en x
for j = 1:M-1
    for i = 1:M-1
        pos = (j-1)*(M-1)+i;                % fila del sistema para el nodo (i,j)
        b(pos) = h^2*f(i*h,j*h);            % lado derecho, ya multiplicado por h^2
        UexMat(i,j) = uex(i*h,j*h);             % referencia para medir el error

        % la fila pos del sistema es la ecuacion del nodo (i,j); todas sus
        % entradas se escriben como (fila pos, columna del vecino)

        % termino diagonal: el propio nodo (i,j)
        k = k+1;  ii(k) = pos;  jj(k) = pos;          ss(k) =  4;

        % vecino izquierdo (i-1,j); si i==1 cae en la frontera x=0
        if(i~=1)
            k = k+1;  ii(k) = pos;  jj(k) = pos-1;      ss(k) = -1;
        end

        % vecino derecho (i+1,j); si i==M-1 cae en la frontera x=1
        if(i~=M-1)
            k = k+1;  ii(k) = pos;  jj(k) = pos+1;      ss(k) = -1;
        end

        % vecino inferior (i,j-1); si j==1 cae en la frontera y=0
        if(j~=1)
            k = k+1;  ii(k) = pos;  jj(k) = pos-(M-1);  ss(k) = -1;
        end

        % vecino superior (i,j+1); si j==M-1 cae en la frontera y=1
        if(j~=M-1)
            k = k+1;  ii(k) = pos;  jj(k) = pos+(M-1);  ss(k) = -1;
        end
    end
end

% control: si el conteo de arriba fuera incorrecto, quedarian ceros al final de
% los vectores (o el indice se habria salido del rango). Conviene verificarlo
assert(k == nnzA, 'se escribieron %d entradas y se habian reservado %d', k, nnzA)

A = sparse(ii,jj,ss,dimA,dimA);             % arma la matriz de una sola vez
uh = A\b;                                   % resuelve el sistema

% visualización
xx = linspace(0,1,M+1);                     % malla completa, con frontera
[X,Y] = meshgrid(xx,xx);
Z = zeros(M+1,M+1);                         % los ceros del borde son la cond. de frontera
% la transpuesta es necesaria porque meshgrid indexa Z(fila,columna) = Z(y,x),
% al reves que UexMat(i,j) = u(x_i,y_j)
Z(2:end-1,2:end-1) = reshape(uh,M-1,M-1)';
surf(X,Y,Z)
xlabel('x')
ylabel('y')

% error en norma infinito
norm(uh-UexMat(:),inf)
