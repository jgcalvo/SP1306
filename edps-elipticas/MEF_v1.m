clear; clc

n = 4;                                  % subdivisiones por lado
h = 1/n;                                 % paso de la malla
f   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
uex = @(x,y) sin(pi*x).*sin(pi*y);

%% ------------------------------------------------------------------
%  1. Malla del cuadrado unitario (guía 1)
[node,elem] = squaremesh_basico(n);
NV = size(node,1);
NT = size(elem,1);
dibujarMalla(node,elem)
etiquetarMalla(node,elem)
[bdNode,~] = nodosFrontera(elem);
freeNode = setdiff(1:NV,bdNode);
%% ------------------------------------------------------------------
%  2. Ensamblar matriz
%     Aloc(i,j) = int_T grad(phi_i) . grad(phi_j) = |T| * (g_i . g_j),
%     bloc(i)   = int_T f phi_i                 ~ |T| f(baricentro)/3,
A = zeros(NV,NV);                        % matriz de rigidez (stiffness)
b = zeros(NV,1);                         % lado derecho
for t = 1:NT
    % --- los tres nodos del elemento y sus coordenadas ---------------
    idx = elem(t,:);                     % indices globales [na nb nc]
    x   = node(idx,1);
    y   = node(idx,2);
    % --- area con signo (positiva porque la malla es antihoraria) ----
    area = 0.5*( (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) );
    % --- gradientes de las tres funciones de forma ------------------
    g = zeros(3,2);
    g(1,:) = [ y(2)-y(3), x(3)-x(2) ] / (2*area);
    g(2,:) = [ y(3)-y(1), x(1)-x(3) ] / (2*area);
    g(3,:) = [ y(1)-y(2), x(2)-x(1) ] / (2*area);
    % --- matrices y vector locales ----------------------------------
    Aloc = zeros(3,3);
    bloc = zeros(3,1);
    xb = (x(1)+x(2)+x(3))/3;             % baricentro del elemento
    yb = (y(1)+y(2)+y(3))/3;
    for i = 1:3
        for j = 1:3
            Aloc(i,j) = area * (g(i,1)*g(j,1) + g(i,2)*g(j,2));
        end
        bloc(i) = area * f(xb,yb) / 3;
    end
    % --- se suma en la matriz global --------------------------------
    for i = 1:3
        for j = 1:3
            A(idx(i),idx(j)) = A(idx(i),idx(j)) + Aloc(i,j);
        end
        b(idx(i)) = b(idx(i)) + bloc(i);
    end
end
% ------------------------------------------------------------------
%  3. Solucion del sistema
u = zeros(NV,1);
u(freeNode) = A(freeNode,freeNode) \ b(freeNode);
% ------------------------------------------------------------------
%  4. Verificacion y grafica
ue  = uex(node(:,1),node(:,2));
err = norm(ue-u,inf)

figure
trisurf(elem, node(:,1), node(:,2), u, 'EdgeColor',[.3 .3 .3]);
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

function dibujarMalla(node,elem)
%DIBUJARMALLA Dibuja una triangulacion en 2D.
clf
patch('Faces',elem,'Vertices',node, ...
    'FaceColor',[0.92 0.94 0.98],'EdgeColor',[0.25 0.25 0.25]);
axis equal tight; box on
end

function etiquetarMalla(node,elem)
%ETIQUETARMALLA Rotula los nodos (azul) y los elementos (rojo).
hold on
plot(node(:,1),node(:,2),'k.','MarkerSize',8);
d = 0.02*max(max(node)-min(node));         % desplazamiento del rotulo
for i = 1:size(node,1)
    text(node(i,1)+d, node(i,2)+d, sprintf('%d',i), ...
        'Color','b','FontSize',9);
end
bc = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
for k = 1:size(elem,1)
    text(bc(k,1), bc(k,2), sprintf('%d',k), 'Color','r', ...
        'FontSize',9,'HorizontalAlignment','center');
end
hold off
end

function [bdNode,bdEdge] = nodosFrontera(elem)
totalEdge  = sort([elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])], 2);
[edge,~,j] = unique(totalEdge,'rows');
conteo     = accumarray(j,1);
bdEdge = edge(conteo==1,:);
bdNode = unique(bdEdge(:));
end
