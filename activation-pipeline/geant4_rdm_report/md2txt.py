"""Convert the bug-report Markdown files to plain ASCII text for Bugzilla.

Usage: python3 md2txt.py */report_*.md   (writes report_*.txt next to each file)
"""
import re, sys, textwrap

UNI = [("ν̄ₑ", "anti_nu_e"), ("νₑ", "nu_e"), ("β⁻", "beta-"), ("β⁺", "beta+"), ("e⁻", "e-"),
       ("e⁺", "e+"), ("β", "beta"), ("α", "alpha"), ("τ", "tau"), ("T½", "T1/2"),
       ("≤", "<="), ("≥", ">="), ("±", "+-"), ("×", "x"), ("−", "-"), ("–", "-"), ("—", "--"),
       ("…", "..."), ("“", '"'), ("”", '"'), ("’", "'"), ("Q_beta-", "Q(beta-)"), ("¹", "(1)")]
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def asc(s):
    for a, b in UNI:
        s = s.replace(a, b)
    s = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+", lambda m: "^" + m.group(0).translate(SUP), s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\1", s)
    s = s.replace("`", "")
    return s


def table(rows):
    cells = [[asc(c.strip()) for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [c for c in cells if not all(re.fullmatch(r":?-+:?", x) for x in c)]
    w = [max(len(r[i]) for r in cells) for i in range(len(cells[0]))]
    out = ["    " + "  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip() for r in cells]
    out.insert(1, "    " + "  ".join("-" * x for x in w))
    return out


def convert(md):
    out, lines, i = [], md.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("<!--"):
            out += ["*** " + asc(ln[4:].rstrip("->").strip()) + " ***"]
        elif ln.startswith("```"):
            i += 1
            while not lines[i].startswith("```"):
                out.append(("    " + asc(lines[i])).rstrip())
                i += 1
        elif ln.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(lines[i]); i += 1
            out += table(rows)
            continue
        elif ln.startswith("# "):
            t = asc(ln[2:])
            out += textwrap.wrap(t, 78)
            out.append("=" * min(78, len(t)))
        elif ln.startswith("## "):
            t = asc(ln[3:]); out += [t, "-" * len(t)]
        elif m := re.match(r"(\s*)([-*]|\d+\.) (.*)", ln):
            ind = m.group(1) + m.group(2) + " "
            out += textwrap.wrap(asc(m.group(3)), 78, initial_indent=ind,
                                 subsequent_indent=" " * len(ind), break_on_hyphens=False)
        elif ln.strip():
            out += textwrap.wrap(asc(ln), 78, break_on_hyphens=False)
        else:
            out.append("")
        i += 1
    txt = "\n".join(out).rstrip() + "\n"
    bad = sorted({c for c in txt if ord(c) > 127})
    if bad:
        sys.exit(f"non-ASCII characters left: {bad}")
    return txt


for p in sys.argv[1:]:
    open(p[:-3] + ".txt", "w").write(convert(open(p).read()))
