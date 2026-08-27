% MDF v2 - Estudio de convergencia, con memoria preasignada
%
% Mismo estudio que MDF_v1_converg.m, con el ensamblaje de MDF_v2_test.m: los
% vectores ii, jj, ss se reservan de entrada con su tamano exacto
%
%     nnz = 5*(M-1)^2 - 4*(M-1)
%
% y se llenan con un contador, en vez de crecer de a un elemento por vez.
%
% Los resultados numericos son identicos a los de la v1 -- misma matriz, misma
% solucion, mismo error. Lo que cambia es el tiempo: al comparar los tic/toc de
% ambas versiones se ve que la v1 se degrada al refinar la malla, porque hacer
% crecer un arreglo dentro de un ciclo cuesta O(N^2), mientras que la v2 se
% mantiene proporcional al tamano del problema.
%
% El error sigue bajando como h^2: cada vez que M se duplica, el error se divide
% entre 4 aproximadamente.

% entradas
uex = @(x,y) sin(pi*x).*sin(pi*y);        % sol. exacta
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y); % lado derecho

err = [];                                 % error de cada malla
hh = [];                                  % paso h de cada malla

for M = 2.^(2:7)                                     % cantidad intervalos
    tic
    h = 1/M;                                         % paso de la malla
    dimA = (M-1)^2;                                  % cantidad de nodos interiores

    % espacio exacto para las tres listas del formato triplete
    nnzA = 5*(M-1)^2 - 4*(M-1);
    ii = zeros(nnzA,1);                              % filas
    jj = zeros(nnzA,1);                              % columnas
    ss = zeros(nnzA,1);                              % valores
    k  = 0;                                          % entradas escritas

    b  = zeros(dimA,1);                              % lado derecho del sistema
    UexMat = zeros(M-1,M-1);                             % sol. exacta en la malla

    % recorrido de los nodos interiores: j avanza en y, i avanza en x
    for j = 1:M-1
        for i = 1:M-1
            pos = (j-1)*(M-1)+i;         % fila del sistema para el nodo (i,j)
            b(pos) = h^2*f(i*h,j*h);     % lado derecho, ya multiplicado por h^2
            UexMat(i,j) = uex(i*h,j*h);      % referencia para medir el error

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
    assert(k == nnzA, 'se escribieron %d entradas y se habian reservado %d', k, nnzA)

    A = sparse(ii,jj,ss,dimA,dimA);      % arma la matriz de una sola vez
    uh = A\b;                            % resuelve el sistema
    toc

    % UexMat(:) apila la matriz por columnas, que es justo el orden
    % pos = (j-1)*(M-1)+i con que se armo uh: ambos vectores son comparables
    hh = [hh, h];
    err = [err, norm(uh-UexMat(:),inf)];
end

%%
% en log-log el error debe quedar paralelo a h^2 (pendiente 2), no a h
loglog(hh,err,'r'), hold on
loglog(hh,hh.^2,'k')
loglog(hh,hh,'b')
legend('err','2','1')
