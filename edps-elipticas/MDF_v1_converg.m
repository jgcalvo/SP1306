% MDF - Estudio de convergencia del metodo de diferencias finitas en 2D
%
% Mismo problema y mismo ensamblaje que MDF_v1_test.m (ver ahi la explicacion
% detallada de la formula de 5 puntos y de la numeracion de los nodos):
%
%     -Lap(u) = f    en  Omega = (0,1)x(0,1)
%           u = 0    sobre la frontera de Omega
%
% Aqui se resuelve el problema para una sucesion de mallas cada vez mas finas,
% M = 4, 8, 16, ..., 128, y se mide el error contra la solucion exacta. El
% objetivo es verificar experimentalmente el orden del metodo.
%
% La teoria dice que la formula de 5 puntos es de orden 2, es decir
% error ~ C*h^2. En escala log-log eso es una recta de pendiente 2, y cada vez
% que h se parte a la mitad el error deberia dividirse entre 4. Por eso al
% final se grafican tambien las rectas h^2 y h: la curva del error debe quedar
% paralela a la primera, no a la segunda.
%
% El tic/toc de cada iteracion muestra el costo de esta version: como ii, jj y
% ss crecen de a un elemento por vez, MATLAB tiene que reservar memoria nueva y
% copiar todo el arreglo en cada paso. El tiempo crece mucho mas rapido que el
% tamano del problema. MDF_v2_converg.m corrige exactamente eso.

% entradas
uex = @(x,y) sin(pi*x).*sin(pi*y);        % sol. exacta
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y); % lado derecho

err = [];                                 % error de cada malla
hh = [];                                  % paso h de cada malla

for M = 2.^(2:7)                                     % cantidad intervalos
    tic
    h = 1/M;                                         % paso de la malla
    dimA = (M-1)^2;                                  % cantidad de nodos interiores
    % vectores para generar A
    ii = []; jj = []; ss = [];
    b = zeros(dimA,1);                               % lado derecho del sistema
    UexMat = zeros(M-1,M-1);                             % sol. exacta en la malla

    % recorrido de los nodos interiores: j avanza en y, i avanza en x
    for j = 1:M-1
        for i = 1:M-1
            pos = (j-1)*(M-1)+i;         % fila del sistema para el nodo (i,j)
            b(pos) = h^2*f(i*h,j*h);     % lado derecho, ya multiplicado por h^2
            UexMat(i,j) = uex(i*h,j*h);      % referencia para medir el error

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
                % i==1: el vecino cae sobre la frontera x=0, donde u=0, asi que
                % no aporta al sistema
                % agregue condicion frontera
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
    A = sparse(ii,jj,ss,dimA,dimA);      % arma la matriz de una sola vez
    %Afull = full(A);
    uh = A\b;                            % resuelve el sistema
    toc
    % % visualización
    % xx = linspace(0,1,M+1);
    % [X,Y] = meshgrid(xx,xx);
    % Z = zeros(M+1,M+1);
    % Z(2:end-1,2:end-1) = reshape(uh,M-1,M-1)';
    % surf(X,Y,Z)
    % xlabel('x')
    % ylabel('y')

    % error norma inf
    % restar uex - uh
    % u1(:) apila la matriz por columnas, que es justo el orden pos = (j-1)*(M-1)+i
    % con que se armo uh, de modo que ambos vectores son comparables entrada a entrada
    hh = [hh, h];
    err = [err, norm(uh-UexMat(:),inf)];
end

%%
% en log-log el error debe quedar paralelo a h^2 (pendiente 2), no a h
loglog(hh,err,'r'), hold on
loglog(hh,hh.^2,'k')
loglog(hh,hh,'b')
legend('err','2','1')
