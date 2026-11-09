# SP1306

Material del curso. Cada tema tiene los scripts originales de MATLAB (`.m`) y su
traducción a notebooks de Python, que se abren en Google Colab con un clic.

## EDPs de orden uno

| Tema | MATLAB | Colab |
|---|---|---|
| Ecuación de transporte, 1D | [`.m`](edps-orden-uno/cap1_transporte1D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_transporte1D.ipynb) |
| Ecuación de transporte, 2D | [`.m`](edps-orden-uno/cap1_transporte2D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_transporte2D.ipynb) |
| Ecuación de Burgers (shock) | [`.m`](edps-orden-uno/cap1_Burgers.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_Burgers.ipynb) |

## EDPs elípticas

Diferencias finitas para $-\Delta u = f$ en el cuadrado unitario, con condiciones
de frontera homogéneas. Las versiones de MATLAB van de lo didáctico a lo óptimo:
la `v1` arma los vectores `ii, jj, ss` por concatenación, la `v2` preasigna el
espacio exacto, y la `v3` construye la matriz sin ciclos, como producto de
Kronecker del laplaciano 1D.

| Tema | MATLAB | Colab |
|---|---|---|
| MDF: solución y visualización | [`v1`](edps-elipticas/MDF_v1_test.m) · [`v2`](edps-elipticas/MDF_v2_test.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-elipticas/MDF_test.ipynb) |
| MDF: estudio de convergencia | [`v1`](edps-elipticas/MDF_v1_converg.m) · [`v2`](edps-elipticas/MDF_v2_converg.m) · [`v3`](edps-elipticas/MDF_v3_converg.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-elipticas/MDF_converg.ipynb) |

## Método de elemento finito

Elementos finitos lineales (P1) para el mismo problema, sobre una triangulación del
cuadrado. La `v1` recorre los elementos uno por uno y usa una matriz llena; la `v2`
arma la matriz rala de forma vectorizada, siguiendo `assembling.m` de las notas.

| Tema | MATLAB | Colab |
|---|---|---|
| MEF P1: solución y convergencia | [`v1`](edps-elipticas/MEF_v1.m) · [`v2`](edps-elipticas/MEF_v2.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-elipticas/MEF.ipynb) |

## EDPs parabólicas

Ecuación del calor con elementos P1 en el espacio. Aquí aparecen las dos matrices
del ensamblaje: la de masa, del término temporal, y la de rigidez, del espacial.
Para avanzar en el tiempo se usa el método $\theta$, que con $\theta = 0, 1, 1/2$
da Euler hacia adelante, Euler hacia atrás y Crank–Nicolson.

| Tema | MATLAB | Colab |
|---|---|---|
| Calor 1D: Euler implícito y explícito | [`.m`](edps-parabolicas/MEF_calor1D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-parabolicas/MEF_calor1D.ipynb) |
| Calor 2D: método $\theta$ y errores en $\Delta t$, $h$ | [`.m`](edps-parabolicas/MEF_calor2D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-parabolicas/MEF_calor2D.ipynb) |

## EDPs hiperbólicas

Ecuación de onda en una dimensión, resuelta en forma cerrada: d'Alembert en toda
la recta, y el método de reflexión cuando el dominio espacial se restringe. Lo
que se ve en las animaciones es qué le pasa a un pulso al llegar al extremo, que
depende solo del coeficiente de reflexión $r = (1-\lambda c)/(1+\lambda c)$:
Neumann lo devuelve sin invertir ($r = 1$), Dirichlet lo invierte ($r = -1$), y
la condición absorbente $\lambda = 1/c$ lo deja salir sin reflejar ($r = 0$).

En dos dimensiones aparece además el elemento finito: la semidiscretización en espacio deja
el sistema $M\bm{U}''+A\bm{U}=\bm{F}$, y lo que queda por decidir es cómo avanzar en el
tiempo. Los tres métodos usuales --- diferencias centradas, Crank--Nicolson y Newmark-$\beta$
--- son en realidad tres puntos de una misma familia de dos parámetros.

Los scripts de 1D escriben sus salidas (PDF o MP4) en el directorio desde el que se corren.

| Tema | MATLAB | Colab |
|---|---|---|
| Onda 2D: $P_1$ en espacio, tres métodos en el tiempo y la CFL | [`.m`](edps-hiperbolicas/demo_onda2d.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-hiperbolicas/demo_onda2d.ipynb) |
| Figuras: reflexión, periodicidad y cuerda pulsada | [`.m`](edps-hiperbolicas/onda1d_figuras.m) | |
| Video: reflexión en $x=0$, Dirichlet contra Neumann | [`.m`](edps-hiperbolicas/onda1d_video.m) | |
| Video: condición de Robin y frontera absorbente | [`.m`](edps-hiperbolicas/onda1d_video_robin.m) | |
| Video: la cuerda pulsada y sus modos | [`.m`](edps-hiperbolicas/onda1d_video_cuerda.m) | |

---
