#!/usr/bin/env python3
"""Inserta (o actualiza) el badge 'Open in Colab' al inicio de cada notebook.

El repo se deduce del remote `origin` de git, asi que no hay nada que editar a
mano cuando se agreguen notebooks nuevos:

    python3 tools/colab_badges.py

Tambien acepta el repo explicito, por si todavia no hay remote:

    python3 tools/colab_badges.py --repo usuario/SP1306 --branch main
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

BADGE_SVG = "https://colab.research.google.com/assets/colab-badge.svg"
RAIZ = pathlib.Path(__file__).resolve().parent.parent


def repo_desde_git():
    """Devuelve 'usuario/repo' leyendo el remote origin, o None si no hay."""
    try:
        url = subprocess.check_output(
            ["git", "-C", str(RAIZ), "remote", "get-url", "origin"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    # soporta git@github.com:usuario/repo.git y https://github.com/usuario/repo.git
    m = re.search(r"github\.com[:/](?P<slug>[^/]+/[^/]+?)(?:\.git)?$", url)
    return m.group("slug") if m else None


def rama_actual():
    try:
        return subprocess.check_output(
            ["git", "-C", str(RAIZ), "rev-parse", "--abbrev-ref", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "main"


def celda_badge(repo, rama, ruta_rel):
    url = f"https://colab.research.google.com/github/{repo}/blob/{rama}/{ruta_rel}"
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [f"[![Open In Colab]({BADGE_SVG})]({url})"],
    }


def es_badge(celda):
    return (celda.get("cell_type") == "markdown"
            and BADGE_SVG in "".join(celda.get("source", [])))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo", help="usuario/repo (por defecto: el remote origin)")
    p.add_argument("--branch", help="rama (por defecto: la rama actual)")
    args = p.parse_args()

    repo = args.repo or repo_desde_git()
    if not repo:
        sys.exit("No hay remote 'origin' todavia. Usá: --repo usuario/SP1306")
    rama = args.branch or rama_actual()

    for nb_path in sorted(RAIZ.rglob("*.ipynb")):
        if ".ipynb_checkpoints" in nb_path.parts:
            continue
        doc = json.loads(nb_path.read_text())
        # la URL de Colab necesita la ruta relativa a la raiz del repo
        ruta_rel = nb_path.relative_to(RAIZ).as_posix()
        nueva = celda_badge(repo, rama, ruta_rel)
        celdas = doc["cells"]

        if celdas and es_badge(celdas[0]):
            celdas[0] = nueva          # actualiza el que ya estaba
            accion = "actualizado"
        else:
            celdas.insert(0, nueva)    # lo agrega arriba de todo
            accion = "agregado"

        nb_path.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
        print(f"badge {accion}: {ruta_rel}")


if __name__ == "__main__":
    main()
