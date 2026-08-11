#!/usr/bin/env python3
"""Ejecuta todos los notebooks de punta a punta y falla si alguno rompe.

Sirve para verificar antes de subir, y es lo que corre GitHub Actions en cada push:

    python3 tools/run_notebooks.py

Las celdas con el tag `solo-colab` (por ejemplo `files.download`) se saltan,
porque dependen del entorno de Colab.
"""
import pathlib
import sys

import nbformat
from nbclient import NotebookClient

RAIZ = pathlib.Path(__file__).resolve().parent.parent
TAG_OMITIR = "solo-colab"


def main():
    notebooks = [p for p in sorted(RAIZ.rglob("*.ipynb"))
                 if ".ipynb_checkpoints" not in p.parts]
    if not notebooks:
        sys.exit("No se encontraron notebooks.")

    fallos = []
    for nb_path in notebooks:
        print(f"\n=== {nb_path.relative_to(RAIZ)}")
        nb = nbformat.read(nb_path, as_version=4)

        omitidas = [c for c in nb.cells
                    if TAG_OMITIR in c.get("metadata", {}).get("tags", [])]
        nb.cells = [c for c in nb.cells if c not in omitidas]
        if omitidas:
            print(f"    ({len(omitidas)} celda(s) '{TAG_OMITIR}' omitidas)")

        # se ejecuta con el cwd en la carpeta del notebook, para que los archivos
        # que genere (videos, imagenes) queden al lado de el y no en la raiz
        cliente = NotebookClient(nb, timeout=900, kernel_name="python3",
                                 resources={"metadata": {"path": str(nb_path.parent)}})
        try:
            cliente.execute()
            print("    OK")
        except Exception as e:
            print(f"    FALLO -> {type(e).__name__}: {e}")
            fallos.append(str(nb_path.relative_to(RAIZ)))

    print("\n" + "=" * 40)
    if fallos:
        print("Fallaron:", ", ".join(fallos))
        sys.exit(1)
    print(f"Los {len(notebooks)} notebooks corrieron sin errores.")


if __name__ == "__main__":
    main()
