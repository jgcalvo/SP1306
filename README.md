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

---
