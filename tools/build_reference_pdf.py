"""Build a single PDF containing all four project markdown documents, verbatim.

Run:  uv run --with markdown --with pypdf python build_reference_pdf.py

Why a script rather than a one-off: this is a point-in-time reference snapshot,
and the repo will change when the D30 -> D52 model lands. Regenerating should be
one command, and the verification step should run every time.

Two markdown/rendering traps this script exists to avoid, both found by the
verification pass rather than by eye:

  * The `attr_list` extension eats `{...}` as an attribute list. The
    hyperparameter grids are written `alpha in {0.01, 0.1, 1, 10, 100}`, so
    enabling it silently deleted every grid from the table -- rendering
    "alpha in" with nothing after it. Do not add attr_list back.
  * `text-transform: uppercase` on table headers rewrites Greek letters:
    eta -> Eta, rho -> Rho. A header reading "pool eta^2" came out as
    "pool Eta^2", which changes what the symbol means.
"""

import re
import subprocess
import unicodedata
from datetime import date
from pathlib import Path

import markdown
from pypdf import PdfReader

ROOT = Path(__file__).parent.parent
OUT_DIR = ROOT / "build"
PDF = OUT_DIR / "pluricon_prototype_reference.pdf"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

PARTS = [
    ("Part I", "Project overview", "README.md",
     "Root README: the question, the data, the pipeline in brief, and headline results."),
    ("Part II", "Findings and literature assessment", "FINDINGS.md",
     "What was found, and which findings recapitulate Jerber et al. (2021) versus add something new."),
    ("Part III", "Modeling methods and full results", "modeling/README.md",
     "The complete method record: CV design, metrics, feature importance, every correction "
     "found by audit, and D30 readiness."),
    ("Part IV", "EDA outputs and provenance", "metadata_eda/README.md",
     "What each EDA script produced, why, and the caveats found along the way."),
]

# NOTE: attr_list is deliberately absent -- see module docstring.
MD_EXTENSIONS = ["tables", "fenced_code", "sane_lists"]

CSS = """
@page { size: A4; margin: 16mm 14mm; }
* { box-sizing: border-box; }
body { font-family:"Charter",Georgia,serif; font-size:9.4pt; line-height:1.5; color:#16181d; margin:0; }
h1,h2,h3,h4 { font-family:"Helvetica Neue",Helvetica,Arial,sans-serif; color:#0d0f14;
  page-break-after:avoid; break-after:avoid; line-height:1.2; }
h1 { font-size:17pt; margin:0 0 .5em; letter-spacing:-.01em; }
h2 { font-size:13pt; margin:1.5em 0 .45em; padding-bottom:.18em; border-bottom:1px solid #d4d7de; }
h3 { font-size:10.6pt; margin:1.15em 0 .35em; }
h4 { font-size:9.6pt; margin:1em 0 .3em; }
p { margin:0 0 .55em; }
ul,ol { margin:.3em 0 .7em; padding-left:1.25em; }
li { margin-bottom:.2em; }
a { color:#1a4f8a; text-decoration:none; }
code { font-family:"SF Mono",Menlo,Consolas,monospace; font-size:.86em; background:#f2f3f6;
  padding:.06em .28em; border-radius:2px; word-break:break-word; }
pre { background:#f7f8fa; border:1px solid #e0e3ea; border-radius:3px; padding:.6em .75em;
  margin:.5em 0 .9em; page-break-inside:avoid; break-inside:avoid; }
pre code { background:none; padding:0; font-size:7.6pt; line-height:1.45;
  white-space:pre-wrap; word-break:break-word; }
blockquote { margin:.6em 0; padding:.45em .8em; border-left:3px solid #9aa3b4; background:#f5f6f9; }
blockquote p:last-child { margin-bottom:0; }
table { border-collapse:collapse; width:100%; margin:.5em 0 1em; font-size:7.2pt;
  page-break-inside:avoid; break-inside:avoid; }
th,td { border:1px solid #d4d7de; padding:3px 5px; text-align:left; vertical-align:top;
  word-break:break-word; }
/* no text-transform here: it rewrites Greek letters in headers (see docstring) */
th { background:#eef0f4; font-family:"Helvetica Neue",Helvetica,sans-serif; font-size:7pt;
  font-weight:600; letter-spacing:.02em; }
tr { page-break-inside:avoid; }
hr { border:none; border-top:1px solid #d4d7de; margin:1.2em 0; }
img { max-width:100%; }
.titlepage { height:245mm; display:flex; flex-direction:column; justify-content:center;
  page-break-after:always; }
.titlepage .kicker { font-family:"Helvetica Neue",sans-serif; font-size:8pt; letter-spacing:.22em;
  text-transform:uppercase; color:#6b7280; margin:0 0 1.4em; }
.titlepage h1 { font-size:30pt; line-height:1.08; margin:0 0 .35em; border:none; }
.titlepage .sub { font-size:12.5pt; color:#4b5261; margin:0 0 2.6em; max-width:36em; }
.titlepage table { width:auto; font-size:9pt; border:none; }
.titlepage td { border:none; padding:3px 26px 3px 0; }
.titlepage td:first-child { font-family:"Helvetica Neue",sans-serif; font-size:7.6pt;
  text-transform:uppercase; letter-spacing:.06em; color:#6b7280; }
.titlepage .note { margin-top:2.4em; font-size:8.4pt; color:#6b7280; border-left:3px solid #c8cdd8;
  padding-left:.9em; max-width:34em; }
.toc { page-break-after:always; }
.toc h1 { border:none; margin-bottom:1.1em; }
.tocpart { margin-bottom:1.5em; page-break-inside:avoid; }
.tocpartname { font-family:"Helvetica Neue",sans-serif; font-size:11pt; font-weight:600; margin:0 0 .1em; }
.tocsrc { font-size:8pt; color:#6b7280; margin:0 0 .45em; }
.toc ul { list-style:none; padding-left:0; margin:0; column-count:2; column-gap:2em; }
.toc li { font-size:8.2pt; color:#3a4150; margin-bottom:.12em; break-inside:avoid; }
.toc li.l3 { padding-left:1.1em; color:#6b7280; font-size:7.8pt; }
.partbreak { page-break-before:always; height:0; }
.parthead { border-bottom:2px solid #16181d; padding-bottom:.7em; margin-bottom:1.4em; }
.partnum { font-family:"Helvetica Neue",sans-serif; font-size:8pt; letter-spacing:.2em;
  text-transform:uppercase; color:#6b7280; margin:0 0 .25em; }
.parttitle { font-size:21pt; margin:0 0 .25em; border:none; }
.partsrc { font-size:8.4pt; color:#6b7280; margin:0 0 .4em; }
.partblurb { font-size:9pt; color:#4b5261; margin:0; max-width:38em; }
.parthead + h1 { font-size:13pt; margin-top:.2em; }
"""

LIGATURES = {"ﬀ": "ff", "ﬁ": "fi", "ﬂ": "fl", "ﬃ": "ffi", "ﬄ": "ffl"}


def alnum_stream(text: str) -> str:
    """Lowercase alphanumeric stream -- immune to rewrapping, hyphenation and ligatures.

    PDF extraction rewraps lines, splits words at hyphens ("cross- validation") and
    emits ligatures, so comparing on raw lines produces false alarms. Reducing both
    sides to bare alphanumerics compares content and nothing else.
    """
    for lig, plain in LIGATURES.items():
        text = text.replace(lig, plain)
    return re.sub(r"[^a-z0-9]", "", unicodedata.normalize("NFKD", text).lower())


def sections(src: str):
    out = []
    for line in src.split("\n"):
        m = re.match(r"^(#{2,3})\s+(.*)$", line)
        if m:
            out.append((len(m.group(1)), re.sub(r"[`*]", "", m.group(2)).strip()))
    return out


def build() -> str:
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    full = subprocess.run(["git", "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()

    body, toc = [], ['<div class="toc"><h1>Contents</h1>']
    for pnum, ptitle, path, blurb in PARTS:
        src = (ROOT / path).read_text()
        slug = re.sub(r"[^a-z0-9]+", "-", ptitle.lower()).strip("-")
        toc.append(f'<div class="tocpart"><p class="tocpartname">{pnum} &middot; {ptitle}</p>'
                   f'<p class="tocsrc"><code>{path}</code></p><ul>')
        toc += [f'<li class="l{lvl}">{name}</li>' for lvl, name in sections(src)]
        toc.append("</ul></div>")
        body.append(f'\n\n<div class="partbreak"></div>\n\n<div class="parthead" id="{slug}">\n'
                    f'<p class="partnum">{pnum}</p><h1 class="parttitle">{ptitle}</h1>\n'
                    f'<p class="partsrc">source: <code>{path}</code></p>\n'
                    f'<p class="partblurb">{blurb}</p></div>\n\n')
        body.append(src)  # verbatim
    toc.append("</div>")

    rendered = markdown.markdown("\n".join(body), extensions=MD_EXTENSIONS)
    title = f"""<div class="titlepage">
  <p class="kicker">Reference snapshot &middot; complete documentation</p>
  <h1>Predicting D52 dopaminergic<br>differentiation efficiency from D11</h1>
  <p class="sub">The full written record of the D11 &rarr; D52 prototype: overview, findings,
  modeling methods and EDA provenance, in one document.</p>
  <table>
    <tr><td>Repository</td><td>pluricon-prototype</td></tr>
    <tr><td>Commit</td><td><code>{full}</code></td></tr>
    <tr><td>Generated</td><td>{date.today().isoformat()}</td></tr>
    <tr><td>Contents</td><td>{' &middot; '.join(p[2] for p in PARTS)}</td></tr>
    <tr><td>Dataset</td><td>Jerber et al. 2021, <i>Nat Genet</i> 53:304&ndash;312 &mdash; 138 iPSC lines</td></tr>
  </table>
  <p class="note">Captured at commit <code>{rev}</code>, before the repository is extended to the
  D30 &rarr; D52 model. Every section of all four source documents is reproduced verbatim;
  nothing is summarised or omitted.</p>
</div>"""

    OUT_DIR.mkdir(exist_ok=True)
    html = ("<!doctype html><html><head><meta charset='utf-8'><title>pluricon-prototype reference"
            f"</title><style>{CSS}</style></head><body>{title}{''.join(toc)}{rendered}</body></html>")
    (OUT_DIR / "combined.html").write_text(html)
    return str(OUT_DIR / "combined.html")


def render(html_path: str) -> None:
    subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={PDF}", f"file://{html_path}"],
                   capture_output=True, check=True)


def verify() -> bool:
    """Every 60-character window of every source must appear in the PDF text."""
    pdf = alnum_stream(" ".join(p.extract_text() or "" for p in PdfReader(PDF).pages))
    print(f"\nPDF: {len(PdfReader(PDF).pages)} pages, {len(pdf):,} alphanumeric chars\n")
    ok = True
    for _, _, path, _ in PARTS:
        raw = (ROOT / path).read_text()
        # Drop fence MARKER lines only (```bash, ```python, ```). The language tag
        # is markdown syntax, not document text, and is correctly not rendered --
        # counting it as content produces false "missing" reports. The fenced
        # body itself is kept and still checked.
        raw = "\n".join(l for l in raw.split("\n") if not l.strip().startswith("```"))
        text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", raw)
        s = alnum_stream(text)
        wins = [s[i:i + 60] for i in range(0, max(1, len(s) - 59), 30)]
        miss = [w for w in wins if w not in pdf]
        ok &= not miss
        print(f"  {path:26s} {len(wins):4d} windows  coverage {100 * (1 - len(miss) / len(wins)):6.2f}%")
        for m in miss[:3]:
            print(f"       MISSING: {m}")
    print("\n  " + ("ALL CONTENT VERIFIED PRESENT" if ok else "CONTENT MISSING -- do not ship"))
    return ok


if __name__ == "__main__":
    render(build())
    raise SystemExit(0 if verify() else 1)
