"""Une la cabecera y portada del informe de P2 con body.tex -> P3_report.tex."""
from pathlib import Path

H = Path(__file__).parent
s = (H.parent.parent / "P2/report/P2_report.tex").read_text(encoding="utf-8")
pre = s[: s.index(r"\begin{document}")]
pre = pre.replace(r"\newcommand{\practiceid}{2}", r"\newcommand{\practiceid}{3}")
pre = pre.replace(r"\newcommand{\practicetitle}{Practice 2}", r"\newcommand{\practicetitle}{Practice 3}")
cover = s[s.index(r"\begin{document}") : s.index(r"\section{Exercise 1}")]
(H / "P3_report.tex").write_text(pre + cover + (H / "body.tex").read_text(encoding="utf-8"), encoding="utf-8")
