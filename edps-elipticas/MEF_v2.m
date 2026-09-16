% MEF v2 - Elementos finitos P1, ensamblaje vectorizado
%
% Mismo problema que MEF_v1.m:
%
%     -Lap(u) = f    en  Omega = (0,1)x(0,1)
%           u = 0    sobre la frontera de Omega
%
% discretizado con elementos finitos lineales (P1) sobre una triangulacion.
% Lo que cambia es COMO se arma la matriz: aca se sigue la tecnica de
% assembling.m de las notas del curso.
%
%
% QUE CAMBIA RESPECTO A LA v1
%
% La v1 recorre los elementos uno por uno y, dentro de cada uno, recorre los
% 3x3 pares de funciones de forma sumando en A(idx(i),idx(j)). Eso son 9*NT
% accesos individuales a una matriz DENSA de NV x NV.
%
% Lo denso es un problema en eficiencia. Con n subdivisiones por lado hay
% NV = (n+1)^2 nodos, asi que A ocupa NV^2 numeros:
%
%     n =  32  ->  A pesa   9 MB
%     n = 128  ->  A pesa  2.2 GB
%
% aunque solo unas 7 entradas por fila sean distintas de cero.
%
% Aca se invierte el orden de los ciclos: en vez de recorrer los NT elementos
% por fuera y los 9 pares por dentro, se recorren los 9 pares por fuera y cada
% uno se calcula para TODOS los elementos a la vez, con operaciones de vectores.
% El ciclo queda de 9 iteraciones, sin importar el tamano de la malla. Los
% valores se acumulan en las listas ii, jj, ss y al final sparse() arma la
% matriz rala de una sola vez.
%
%
% LA IDENTIDAD QUE PERMITE VECTORIZAR
%
% Sea T un triangulo de vertices p1, p2, p3 y area |T|, y sea
%
%     v_i = el vector del lado OPUESTO al vertice i
%
% (o sea v_1 = p3-p2, v_2 = p1-p3, v_3 = p2-p1). La funcion de forma P1
% asociada al vertice i tiene gradiente constante, y ese gradiente es el lado
% opuesto rotado 90 grados y dividido entre 2|T|:
%
%     grad(phi_i) = rot90(v_i) / (2|T|)
%
% Como una rotacion no cambia los productos punto,
%
%     grad(phi_i) . grad(phi_j) = (v_i . v_j) / (4|T|^2)
%
% y al integrar sobre T (el integrando es constante) queda
%
%     int_T grad(phi_i).grad(phi_j) = (v_i . v_j) / (4|T|)
%
% Esa es la formula que se evalua para todos los elementos. Notese que
% solo se calculan restas de coordenadas y un producto punto.
%
% El area se calcula con valor absoluto, asi que la formula no depende de la
% orientacion de los triangulos. La v1 usa el area con signo y solo funciona si
% la malla esta orientada en sentido antihorario.
%

clear; clc

n = 32;                                  % subdivisiones por lado
h = 1/n;                                 % paso de la malla
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
uex = @(x,y) sin(pi*x).*sin(pi*y);

%% ------------------------------------------------------------------
%  1. Malla del cuadrado unitario
[node,elem] = squaremesh_basico(n);
NV = size(node,1);
NT = size(elem,1);
[bdNode,~] = nodosFrontera(elem);
freeNode = setdiff(1:NV,bdNode);

%% ------------------------------------------------------------------
%  2. Ensamblar, siguiendo assembling.m
tic
% ve(:,:,i) = lado opuesto al vertice i, para todos los elementos a la vez
ve = zeros(NT,2,3);
ve(:,:,1) = node(elem(:,3),:) - node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:) - node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:) - node(elem(:,1),:);

% area de cada triangulo, a partir del producto cruz de dos lados
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2) + ve(:,2,3).*ve(:,1,2));

% listas del formato triplete: 9 entradas (3x3) por cada elemento
ii = zeros(9*NT,1);  jj = zeros(9*NT,1);  sA = zeros(9*NT,1);

index = 0;
for i = 1:3
    for j = 1:3
        % cada asignacion llena NT posiciones de un solo golpe
        ii(index+1:index+NT) = elem(:,i);        % fila global
        jj(index+1:index+NT) = elem(:,j);        % columna global
        sA(index+1:index+NT) = dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        index = index + NT;
    end
end
% un mismo par (fila,columna) aparece varias veces, una por cada elemento que
% comparte esos dos nodos: sparse SUMA los valores repetidos, que es justo el
% ensamblaje
A = sparse(ii,jj,sA,NV,NV);

% lado derecho: misma regla del baricentro que la v1, tambien vectorizada
xb = (node(elem(:,1),1) + node(elem(:,2),1) + node(elem(:,3),1))/3;
yb = (node(elem(:,1),2) + node(elem(:,2),2) + node(elem(:,3),2))/3;
aporte = area.*f(xb,yb)/3;               % lo que cada elemento da a c/u de sus 3 nodos
b = accumarray(elem(:), repmat(aporte,3,1), [NV 1]);
toc

%% ------------------------------------------------------------------
%  3. Solucion del sistema
u = zeros(NV,1);
u(freeNode) = A(freeNode,freeNode) \ b(freeNode);

%% ------------------------------------------------------------------
%  4. Verificacion y grafica
ue  = uex(node(:,1),node(:,2));
err = norm(ue-u,inf)

figure
trisurf(elem, node(:,1), node(:,2), u, 'EdgeColor','none');
colorbar
xlabel('x'); ylabel('y'); zlabel('u_h');

%% funciones auxiliares
function [node,elem] = squaremesh_basico(M)
%SQUAREMESH_BASICO Malla estructurada de (0,1)^2, M subintervalos por lado.
%   node: (M+1)^2 x 2   elem: 2*M^2 x 3, orientados en sentido antihorario.
[X,Y] = meshgrid(linspace(0,1,M+1));
node  = [X(:), Y(:)];   % el nodo (x_i,y_j) ocupa la fila (i-1)*(M+1)+j
elem  = zeros(2*M^2,3);
k = 0;
for i = 1:M                    % columnas (direccion x)
    for j = 1:M                % filas    (direccion y)
        n1 = (i-1)*(M+1) + j;  % (x_i    , y_j    )
        n2 = n1 + (M+1);       % (x_{i+1}, y_j    )
        n3 = n2 + 1;           % (x_{i+1}, y_{j+1})
        n4 = n1 + 1;           % (x_i    , y_{j+1})
        elem(k+1,:) = [n1 n2 n3];   % antihorario
        elem(k+2,:) = [n1 n3 n4];   % antihorario
        k = k + 2;
    end
end
end

function [bdNode,bdEdge] = nodosFrontera(elem)
%NODOSFRONTERA Nodos y aristas de la frontera de una triangulacion.
%   Una arista es de frontera si pertenece a un solo triangulo.
totalEdge  = sort([elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])], 2);
[edge,~,j] = unique(totalEdge,'rows');
conteo     = accumarray(j,1);
bdEdge = edge(conteo==1,:);
bdNode = unique(bdEdge(:));
end
