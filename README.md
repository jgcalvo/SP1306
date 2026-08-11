# SP1306

Material del curso. Cada tema tiene los scripts originales de MATLAB (`.m`) y su
traducción a notebooks de Python, que se abren en Google Colab con un clic.

## EDPs de orden uno

| Tema | MATLAB | Colab |
|---|---|---|
| Ecuación de transporte, 1D | [`.m`](edps-orden-uno/cap1_transporte1D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_transporte1D.ipynb) |
| Ecuación de transporte, 2D | [`.m`](edps-orden-uno/cap1_transporte2D.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_transporte2D.ipynb) |
| Ecuación de Burgers (shock) | [`.m`](edps-orden-uno/cap1_Burgers.m) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jgcalvo/SP1306/blob/main/edps-orden-uno/cap1_Burgers.ipynb) |

## Método de elemento finito

_Pendiente._

---

## Cómo abrirlos en Colab

Cada notebook lleva arriba un badge **Open in Colab**: al hacerle clic desde GitHub,
Colab lo abre directo, sin descargar nada. Ese es todo el flujo.

Para que los cambios que hagas en Colab vuelvan al repo: en Colab,
`Archivo > Guardar una copia en GitHub` (la primera vez pide autorizar la app de
Colab en tu cuenta). Si preferís no autorizarla, bajá el `.ipynb` y lo subimos desde acá.

No hace falta instalar nada: Colab ya trae `numpy`, `matplotlib` y `ffmpeg`.

## Flujo de trabajo

```
.m nuevo en la carpeta del tema  →  traducción a .ipynb  →  push  →  badge → Colab
```

Después de agregar un notebook nuevo (los scripts recorren las subcarpetas solos):

```bash
python3 tools/colab_badges.py    # pone/actualiza el badge (lee el remote de git)
python3 tools/run_notebooks.py   # verifica que todos corren sin errores
```

`colab_badges.py` deduce el repo del remote `origin`, así que no hay nada que editar
a mano. `run_notebooks.py` salta las celdas con el tag `solo-colab` (las que usan
`google.colab`, que no existe fuera de Colab).

## ¿Correrlos en GitHub?

Tres cosas distintas que conviene no confundir:

- **Ver.** GitHub renderiza los `.ipynb` solo. Si el notebook se guardó con salidas,
  se ven los gráficos sin ejecutar nada. Es lectura, no ejecución.
- **Ejecutar de verdad, interactivo.** [GitHub Codespaces](https://github.com/features/codespaces)
  te levanta VS Code con Jupyter en el navegador (las cuentas gratuitas traen horas
  incluidas por mes). Sirve, pero para animaciones Colab es más cómodo y no consume cuota.
- **Ejecutar automático, sin interacción.** Es lo que hace `.github/workflows/notebooks.yml`:
  en cada push corre los tres notebooks en un runner de Ubuntu y falla si alguno rompe.
  Es un chequeo de que el código sigue funcionando, no una forma de trabajar.

En resumen: **GitHub guarda y versiona, Colab ejecuta.** Los Actions te avisan si algo
se rompió.

## Diferencias respecto a MATLAB

- **Animaciones.** El patrón `plot` + `pause(.1)` dentro de un `for` no funciona en
  Colab (la salida de una celda se dibuja al final). En su lugar se usa
  `matplotlib.animation.FuncAnimation` y se muestra con `HTML(ani.to_jshtml())`,
  que da un reproductor con controles dentro del notebook.
- **Índices.** Python arranca en 0 y `range` excluye el extremo derecho, así que
  `for j = 1:5:numel(tt)` se vuelve `range(0, len(tt), 5)`.
- **`linspace(a,b)`.** En MATLAB son 100 puntos por defecto; en NumPy son 50, por eso
  el `100` va explícito.
- **Malla de la 2D.** El script original usa `1e3 x 1e3` puntos. Dibujar un millón de
  puntos por cuadro en matplotlib es demasiado lento, así que se usan 200 por dirección
  (variable `N`). Se puede subir si querés más detalle.
- **Video.** `VideoWriter('transporte.avi')` se reemplaza por `ani.save('transporte.mp4')`.
  El archivo queda en el disco temporal de Colab y se borra al cerrar la sesión: hay una
  celda con `files.download(...)` para bajarlo.
