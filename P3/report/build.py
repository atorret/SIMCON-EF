"""Cabecera y portada del informe de P2 + \\input{body.tex} -> P3_report.tex, y compila con SyncTeX.

El cuerpo se incluye con \\input (no se copia), asi Ctrl+clic en el PDF salta a body.tex.
Uso:  python build.py
"""
import subprocess
from pathlib import Path

H = Path(__file__).parent
s = (H.parent.parent / "P2/report/P2_report.tex").read_text(encoding="utf-8")
pre = s[: s.index(r"\begin{document}")]
pre = pre.replace(r"\newcommand{\practiceid}{2}", r"\newcommand{\practiceid}{3}")
pre = pre.replace(r"\newcommand{\practicetitle}{Practice 2}", r"\newcommand{\practicetitle}{Practice 3}")
cover = s[s.index(r"\begin{document}") : s.index(r"\section{Exercise 1}")]
(H / "P3_report.tex").write_text(pre + cover + "\\input{body.tex}\n", encoding="utf-8")

subprocess.run(
    ["latexmk", "-pdf", "-synctex=1", "-interaction=nonstopmode", "-halt-on-error", "P3_report.tex"],
    cwd=H, check=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
)
print("PDF:", H / "P3_report.pdf")
