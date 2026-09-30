#!/usr/bin/env python3
"""
defensome - comparative chemical defensome annotation from proteomes.

  defensome.py scan     --proteomes DIR --pfam Pfam-A.hmm --out DIR [--threads N]
  defensome.py annotate --out DIR [--map defensome_map.tsv]
  defensome.py report   --out DIR [--metadata meta.tsv --group-by COLUMN]
  defensome.py all      --proteomes DIR --pfam Pfam-A.hmm --out DIR [...]

Only hard dependencies: HMMER 3 (for `scan`), Python 3.9+, pandas, numpy.
matplotlib is optional; figures are skipped with a warning if it is missing.
"""
import argparse, glob, gzip, math, os, re, shutil, subprocess, sys

# numpy/pandas are needed by annotate and report, but NOT by scan. Importing
# them lazily means `scan` runs in a bare HMMER module environment, which
# matters on clusters where the HMMER and SciPy toolchains cannot coexist.
try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    np = pd = None


def need_pandas():
    if pd is None:
        sys.exit("ERROR: this step needs pandas and numpy.\n"
                 "  module load SciPy-bundle    (or: pip install --user pandas numpy)")

__version__ = "3.1.0"

DOMTBL_COLS = ["target","t_acc","tlen","query","q_acc","qlen","e_full","score_full",
               "bias_full","dom_n","dom_of","c_evalue","i_evalue","score_dom","bias_dom",
               "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc"]

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MAP = os.path.join(_HERE, "defensome_map.tsv")

# The dashboard template, its JavaScript and the default map are embedded below
# so that defensome.py works on its own from any directory. A file of the same
# name sitting next to the script wins, so you can still edit any of them:
#   python3 defensome.py assets --write .
_ASSETS = {
    "defensome_map.tsv": (
        "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"
        "OpA1yaNkcTVaPK+ael/A6yX+u6Bq4/8HvpLgKA=="
    ),
    "dashboard.html": (
        "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"
    ),
    "dashboard.js": (
        "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"
    ),
}


def asset(name):
    """Sibling file if present, otherwise the copy embedded in this script."""
    p = os.path.join(_HERE, name)
    if os.path.exists(p):
        return open(p).read()
    import base64, zlib
    if name not in _ASSETS:
        sys.exit(f"ERROR: no embedded asset '{name}' and no file at {p}")
    return zlib.decompress(base64.b64decode(_ASSETS[name])).decode()


def cmd_assets(a):
    """Write the embedded assets out so they can be edited."""
    d = a.write or "."
    os.makedirs(d, exist_ok=True)
    for name in ("defensome_map.tsv", "dashboard.html", "dashboard.js"):
        dest = os.path.join(d, name)
        if os.path.exists(dest) and not a.force:
            log(f"exists, skipping (use --force): {dest}"); continue
        import base64, zlib
        with open(dest, "w") as fh:
            fh.write(zlib.decompress(base64.b64decode(_ASSETS[name])).decode())
        log(f"wrote {dest}")


# --------------------------------------------------------------------------- #
def log(msg):
    print(f"[defensome] {msg}", file=sys.stderr, flush=True)


def species_list(proteome_dir):
    fa = sorted(glob.glob(os.path.join(proteome_dir, "*.faa")) +
                glob.glob(os.path.join(proteome_dir, "*.fasta")) +
                glob.glob(os.path.join(proteome_dir, "*.fa")))
    if not fa:
        sys.exit(f"ERROR: no .faa/.fa/.fasta files in {proteome_dir}")
    return [(os.path.splitext(os.path.basename(f))[0], f) for f in fa]


def read_map(path):
    need_pandas()
    import io
    src = path if os.path.exists(path) else None
    if src is None and os.path.abspath(path) == os.path.abspath(DEFAULT_MAP):
        df = pd.read_csv(io.StringIO(asset("defensome_map.tsv")), sep="\t", comment="#")
    elif src is None:
        sys.exit(f"ERROR: map file not found: {path}")
    else:
        df = pd.read_csv(src, sep="\t", comment="#")
    need = {"category","family","pfam_ids","rule","min_cov","min_len","tier"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} is missing columns: {sorted(missing)}")
    df["pfam_ids"] = df.pfam_ids.str.replace(" ", "").str.split(",")
    # optional columns (older maps do not have them)
    if "rescue" not in df.columns:
        df["rescue"] = "no"
    df["rescue"] = df.rescue.fillna("no").astype(str).str.strip().str.lower()
    if "validate" not in df.columns:
        df["validate"] = ""
    df["validate"] = df.validate.fillna("").astype(str).str.strip()
    bad = df[~df.rule.isin(["ANY", "ALL"])]
    if len(bad):
        sys.exit(f"ERROR: rule must be ANY or ALL; bad rows: {list(bad.family)}")
    return df


# =========================================================================== #
#  Pfam HMM file utilities (pure Python, no HMMER needed)
# =========================================================================== #
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_VERSION_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz"


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def iter_hmm_records(path):
    """Yield (header_dict, record_text) for each HMM in a HMMER3 text file."""
    buf, hdr = [], {}
    with _open_text(path) as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("//"):
                yield hdr, "".join(buf)
                buf, hdr = [], {}
                continue
            if len(line) > 6 and line[:6].rstrip() in ("NAME", "ACC", "LENG", "GA", "DESC"):
                key = line[:6].rstrip()
                if key not in hdr:
                    hdr[key] = line[6:].strip()


def hmm_index(path):
    """accession (no version) -> dict(name, leng, ga_seq, ga_dom, desc)"""
    idx = {}
    for h, _ in iter_hmm_records(path):
        acc = h.get("ACC", "").split(".")[0]
        if not acc:
            continue
        ga = h.get("GA", "").rstrip(";").split()
        idx[acc] = {"name": h.get("NAME", ""), "leng": int(h.get("LENG", "0") or 0),
                    "ga_seq": float(ga[0]) if len(ga) > 0 else None,
                    "ga_dom": float(ga[1]) if len(ga) > 1 else None,
                    "desc": h.get("DESC", "")}
    return idx


def write_hmm_subset(src, accessions, dest):
    """Copy only the HMMs whose accession is in `accessions`. Atomic: parallel
    array tasks may race to build the same file; identical content, os.replace."""
    want = {a.split(".")[0] for a in accessions}
    found = set()
    d = os.path.dirname(os.path.abspath(dest))
    os.makedirs(d, exist_ok=True)
    tmp = f"{dest}.{os.getpid()}.tmp"
    with open(tmp, "w") as out:
        for h, rec in iter_hmm_records(src):
            acc = h.get("ACC", "").split(".")[0]
            if acc in want:
                out.write(rec)
                found.add(acc)
    os.replace(tmp, dest)
    return found


def _map_rows_plain(mapfile):
    """Minimal map reader for steps that must run without pandas (scan)."""
    import csv, io
    if os.path.exists(mapfile):
        text = open(mapfile).read()
    else:
        text = asset("defensome_map.tsv")
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    rows = list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))
    for r in rows:
        r["pfam_ids"] = [x.strip() for x in (r.get("pfam_ids") or "").split(",") if x.strip()]
        r["rescue"] = (r.get("rescue") or "no").strip().lower()
    return rows


def _download(url, dest):
    import urllib.request
    log(f"downloading {url}")
    tmp = dest + ".part"
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        got, last = 0, 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if total and got - last > total // 20:
                last = got
                print(f"   {got/1e6:8.1f} / {total/1e6:.1f} MB", file=sys.stderr)
    os.replace(tmp, dest)


def cmd_setup(a):
    """Prepare the HMM database: extract the defensome families from Pfam-A.

    Scanning 30,000 proteins against all ~25,000 Pfam families takes hours per
    proteome. The map needs about 45 of them. Family calls are identical either
    way, because hmmsearch scores each HMM independently and --cut_ga uses
    per-family bit-score thresholds, so the subset gives the same counts in a
    fraction of the time. The full Pfam-A is only needed for the optional
    architecture pass, which reports every other domain on defensome proteins.

    Writes <db-dir>/defensome.hmm and <db-dir>/MANIFEST.tsv. The manifest says
    which accessions your Pfam release actually contains: an accession missing
    from it is the first thing to rule out when a family comes back empty."""
    rows = _map_rows_plain(a.map)
    accs = sorted({x for r in rows for x in r["pfam_ids"]})
    os.makedirs(a.db_dir, exist_ok=True)
    src = a.pfam
    if a.download:
        src = os.path.join(a.db_dir, "Pfam-A.hmm.gz")
        if not os.path.exists(src) or a.force:
            _download(a.url or PFAM_URL, src)
        try:
            vgz = os.path.join(a.db_dir, "Pfam.version.gz")
            _download(PFAM_VERSION_URL, vgz)
            with gzip.open(vgz, "rt") as fh:
                open(os.path.join(a.db_dir, "PFAM_VERSION"), "w").write(fh.read())
        except Exception as e:
            log(f"could not fetch Pfam.version ({e}); release not recorded")
    if not src or not os.path.exists(src):
        sys.exit("ERROR: give --pfam /path/to/Pfam-A.hmm[.gz], or --download")

    log(f"indexing {src} (reads the whole file once; a few minutes for full Pfam-A)")
    dest = os.path.join(a.db_dir, "defensome.hmm")
    found = write_hmm_subset(src, accs, dest)
    idx = hmm_index(dest)
    fam_of = {}
    for r in rows:
        for x in r["pfam_ids"]:
            fam_of.setdefault(x, []).append(r["family"])
    man = os.path.join(a.db_dir, "MANIFEST.tsv")
    with open(man, "w") as fh:
        fh.write("pfam_id\tfamilies\tfound\tname\tlength\tga_seq\tga_dom\tshort_domain\tdescription\n")
        for x in accs:
            i = idx.get(x, {})
            short = "yes" if i and i["leng"] and i["leng"] < 60 else ("no" if i else "")
            fh.write(f"{x}\t{','.join(fam_of.get(x, []))}\t{'yes' if x in found else 'NO'}\t"
                     f"{i.get('name','')}\t{i.get('leng','')}\t{i.get('ga_seq','')}\t"
                     f"{i.get('ga_dom','')}\t{short}\t{i.get('desc','')}\n")
    missing = [x for x in accs if x not in found]
    print(f"\n=== {len(found)} of {len(accs)} map accessions found in {os.path.basename(src)} ===")
    for x in accs:
        i = idx.get(x)
        if i:
            flag = "  SHORT DOMAIN" if i["leng"] < 60 else ""
            print(f"  {x}  {i['name']:18s} length {i['leng']:4d}  GA {i['ga_seq']}{flag}")
    if missing:
        print("\nMISSING from this Pfam release (these families cannot be called):")
        for x in missing:
            print(f"  {x}  used by {', '.join(fam_of.get(x, []))}")
        print("Check the accession on https://www.ebi.ac.uk/interpro/ ; Pfam retires")
        print("and merges families between releases.")
    if a.download and a.keep_full:
        full = os.path.join(a.db_dir, "Pfam-A.hmm")
        with gzip.open(src, "rb") as fi, open(full, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        log(f"kept full Pfam-A for the optional architecture pass: {full}")
    ver = os.path.join(a.db_dir, "PFAM_VERSION")
    if os.path.exists(ver):
        print("\n" + open(ver).read().strip())
    print(f"\nwrote {dest}  ({os.path.getsize(dest)/1e6:.1f} MB, {len(found)} HMMs)")
    print(f"      {man}")
    print(f"\nscan with:  python3 defensome.py scan --pfam {dest} ...")


# --------------------------------------------------------------------------- #
def _gz_move(tmp, dest):
    with open(tmp, "rb") as fi, gzip.open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(tmp)


def cmd_scan(a):
    """hmmsearch per proteome, in up to three passes.

    1. Main pass, --cut_ga: Pfam's curated per-family bit-score thresholds.
       An E-value cutoff would depend on proteome size and so would not be
       comparable across the proteomes being compared. -> hmmsearch/
    2. Rescue pass, for families marked rescue=yes in the map: the same HMMs
       with an E-value cutoff instead of --cut_ga. Very short domains have few
       positions to accumulate score, so a genuine but divergent copy can sit
       below a threshold tuned on other taxa. These hits are kept apart and
       must pass the map's validate rules before they count. -> rescue/
    3. Architecture pass, only with --full-pfam: every Pfam domain on the
       proteins that hit a defensome family, for the architecture strings.
       Restricted to those proteins, so it takes minutes. -> hmmsearch_full/
    """
    hs = tool("hmmsearch")
    if not hs:
        sys.exit("ERROR: hmmsearch not found.\n"
                 "  module load HMMER    (or set DEFENSOME_HMMSEARCH=/path/to/hmmsearch)\n"
                 "  python3 defensome.py doctor   shows everything that is missing")
    outdir = os.path.join(a.out, "hmmsearch"); os.makedirs(outdir, exist_ok=True)
    # record which database was searched, so qc can find its MANIFEST later
    with open(os.path.join(a.out, "scan_db.txt"), "w") as fh:
        fh.write(os.path.abspath(a.pfam) + "\n")
    todo = species_list(a.proteomes)
    if a.species:
        todo = [(s_, f) for s_, f in todo if s_ == a.species]
        if not todo:
            sys.exit(f"ERROR: species '{a.species}' not found in {a.proteomes}")

    # rescue HMMs: built once from the scan database, atomically
    rescue_hmm, rescue_accs = None, set()
    if not a.no_rescue:
        rows = _map_rows_plain(a.map)
        rescue_accs = {x for r in rows if r["rescue"] == "yes" for x in r["pfam_ids"]}
        if rescue_accs:
            rdir = os.path.join(a.out, "rescue"); os.makedirs(rdir, exist_ok=True)
            rescue_hmm = os.path.join(rdir, "rescue_families.hmm")
            if not os.path.exists(rescue_hmm) or a.force:
                got = write_hmm_subset(a.pfam, rescue_accs, rescue_hmm)
                lost = sorted(rescue_accs - got)
                if lost:
                    log(f"WARNING: rescue accessions not in {a.pfam}: {lost}")
                if not got:
                    rescue_hmm = None

    for sp, faa in todo:
        dest = os.path.join(outdir, f"{sp}.domtblout.gz")
        if os.path.exists(dest) and not a.force:
            log(f"{sp}: main pass exists, skipping (use --force to redo)")
        else:
            log(f"{sp}: hmmsearch --cut_ga")
            tmp = dest[:-3]
            subprocess.run([hs, "--cut_ga", "--cpu", str(a.threads),
                            "--domtblout", tmp, "-o", os.devnull, a.pfam, faa], check=True)
            _gz_move(tmp, dest)

        if rescue_hmm:
            rdest = os.path.join(a.out, "rescue", f"{sp}.domtblout.gz")
            if not os.path.exists(rdest) or a.force:
                log(f"{sp}: rescue pass, E <= {a.rescue_evalue}, no --cut_ga")
                tmp = rdest[:-3]
                subprocess.run([hs, "-E", str(a.rescue_evalue), "--domE", str(a.rescue_evalue),
                                "--cpu", str(a.threads), "--domtblout", tmp,
                                "-o", os.devnull, rescue_hmm, faa], check=True)
                _gz_move(tmp, rdest)

        if a.full_pfam:
            fdir = os.path.join(a.out, "hmmsearch_full"); os.makedirs(fdir, exist_ok=True)
            fdest = os.path.join(fdir, f"{sp}.domtblout.gz")
            if not os.path.exists(fdest) or a.force:
                cand = set()
                with gzip.open(dest, "rt") as fh:
                    for line in fh:
                        if not line.startswith("#") and line.strip():
                            cand.add(line.split(None, 1)[0])
                sub = os.path.join(fdir, f"{sp}.candidates.faa")
                n = 0
                with open(sub, "w") as fo:
                    for h, seq in read_fasta(faa):
                        if h.split()[0] in cand:
                            fo.write(f">{h}\n{seq}\n"); n += 1
                log(f"{sp}: architecture pass on {n} candidate proteins against full Pfam")
                tmp = fdest[:-3]
                subprocess.run([hs, "--cut_ga", "--cpu", str(a.threads), "--domtblout", tmp,
                                "-o", os.devnull, a.full_pfam, sub], check=True)
                _gz_move(tmp, fdest)
                os.remove(sub)
    log("scan done")


# --------------------------------------------------------------------------- #
def read_domtbl(path):
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 22)
            if len(f) >= 22:
                rows.append(f[:22])
    if not rows:
        return pd.DataFrame(columns=DOMTBL_COLS + ["pfam_id"])
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)
    for c in ("tlen", "qlen", "hmm_from", "hmm_to"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["pfam_id"] = df.q_acc.str.split(".").str[0]
    return df.dropna(subset=["tlen", "qlen"])


def domain_coverage(df):
    """Per (protein, pfam), sum non-overlapping HMM coverage.
    Repeat domains and split hits are merged rather than double counted."""
    out = []
    for (prot, pf), g in df.groupby(["target", "pfam_id"], sort=False):
        iv = sorted(zip(g.hmm_from.astype(int), g.hmm_to.astype(int)))
        merged, cs, ce = [], None, None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s <= ce + 1:
                ce = max(ce, e)
            else:
                merged.append((cs, ce)); cs, ce = s, e
        if cs is not None:
            merged.append((cs, ce))
        qlen = float(g.qlen.iloc[0])
        out.append((prot, pf, int(g.tlen.iloc[0]),
                    sum(e - s + 1 for s, e in merged) / qlen if qlen else 0.0))
    return pd.DataFrame(out, columns=["protein", "pfam_id", "prot_len", "hmm_cov"])


def count_proteins(faa):
    n = 0
    op = gzip.open if faa.endswith(".gz") else open
    with op(faa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


_VAL_RE = re.compile(r"^(len|cys|copies)\s*(<=|>=|<|>|==)\s*([0-9.]+)$")


def _parse_validate(spec, family=""):
    rules = []
    for part in str(spec or "").split(";"):
        part = part.strip()
        if not part:
            continue
        m = _VAL_RE.match(part)
        if not m:
            sys.exit(f"ERROR: bad validate rule '{part}' for {family}. "
                     "Use len<=N, len>=N, cys>=F or copies>=N, separated by ';'")
        rules.append((m.group(1), m.group(2), float(m.group(3))))
    return rules


def _check_rules(rules, vals):
    """vals: dict key -> number or None. Returns (ok, reason)."""
    ops = {"<=": lambda x, y: x <= y, ">=": lambda x, y: x >= y, "<": lambda x, y: x < y,
           ">": lambda x, y: x > y, "==": lambda x, y: x == y}
    for k, op, v in rules:
        x = vals.get(k)
        if x is None:
            return False, f"{k} unavailable (proteome not given)"
        if not ops[op](x, v):
            return False, f"{k}={x:.3g} fails {k}{op}{v:g}"
    return True, "passes"


MT_RULE = {"min_len": 25, "max_len": 120, "min_cys": 8, "min_cys_frac": 0.18,
           "min_motifs": 4, "max_aromatic": 1}


def mt_features(seq):
    u = seq.upper().rstrip("*")
    L = len(u)
    c = u.count("C")
    return {"length": L, "n_cys": c, "cys_fraction": round(c / L, 3) if L else 0.0,
            "n_motifs": len(re.findall(r"(?=C.C|C..C)", u)),
            "n_aromatic": sum(u.count(x) for x in "FWY")}


def _is_mt_like(seq):
    """Composition screen for metallothionein-like proteins.

    No Pfam model covers insect metallothioneins (InterPro IPR000966, family 5,
    has no Pfam member), so they are found by what defines them: short, rich
    in cysteine arranged in C-x-C and C-x-x-C motifs, and devoid of aromatic
    residues. The aromatic rule does the discriminating work. Cysteine alone
    admits defensins, knottins and chitin-binding peptides, which are also
    short and cysteine-rich but carry F, W or Y; an earlier version of this
    screen used cysteine only and accepted such a decoy while rejecting a real
    metallothionein at 23% cysteine. Heuristic: reported apart from HMM calls."""
    f = mt_features(seq)
    r = MT_RULE
    return (r["min_len"] <= f["length"] <= r["max_len"] and f["n_cys"] >= r["min_cys"]
            and f["cys_fraction"] >= r["min_cys_frac"] and f["n_motifs"] >= r["min_motifs"]
            and f["n_aromatic"] <= r["max_aromatic"])


def _copies(df):
    """distinct non-overlapping target-sequence segments hit by one HMM"""
    iv = sorted(zip(pd.to_numeric(df.ali_from).astype(int), pd.to_numeric(df.ali_to).astype(int)))
    n, ce = 0, -1
    for s_, e_ in iv:
        if s_ > ce:
            n += 1
            ce = e_
        else:
            ce = max(ce, e_)
    return n


def cmd_annotate(a):
    need_pandas()
    dm = read_map(a.map)
    scan_dir = os.path.join(a.out, "hmmsearch")
    files = sorted(glob.glob(os.path.join(scan_dir, "*.domtblout*")))
    if not files:
        sys.exit(f"ERROR: no domtblout files in {scan_dir}. Run `scan` first.")
    calls_dir = os.path.join(a.out, "gene_calls"); os.makedirs(calls_dir, exist_ok=True)
    rescue_fams = dm[dm.rescue == "yes"]
    has_mt = "Metallothionein" in set(dm.family)
    faa_of = dict(species_list(a.proteomes)) if a.proteomes and os.path.isdir(a.proteomes) else {}

    prot_totals, stat_rows, comp_rows, resc_rows, all_calls = {}, [], [], [], []
    for path in files:
        sp = re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(path))
        cov = domain_coverage(read_domtbl(path))
        rows = []
        if not cov.empty:
            has = cov.groupby("protein").pfam_id.apply(set).to_dict()
            for _, fam in dm.iterrows():
                sel = cov[(cov.pfam_id.isin(fam.pfam_ids)) &
                          (cov.hmm_cov >= fam.min_cov) & (cov.prot_len >= fam.min_len)]
                for prot, g in sel.groupby("protein"):
                    if fam.rule == "ALL" and not set(fam.pfam_ids).issubset(has.get(prot, set())):
                        continue
                    rows.append((sp, fam.category, fam.family, fam.tier, prot,
                                 int(g.prot_len.iloc[0]), round(float(g.hmm_cov.max()), 3), "GA"))
        else:
            log(f"WARNING: {sp} produced no Pfam hits")
        c = pd.DataFrame(rows, columns=["species", "category", "family", "tier",
                                        "protein", "prot_len", "hmm_cov", "evidence"])
        called = set(zip(c.family, c.protein))

        # rescue candidates from the below-threshold pass
        cands = []
        rpath = os.path.join(a.out, "rescue", f"{sp}.domtblout.gz")
        if len(rescue_fams) and os.path.exists(rpath):
            rd = read_domtbl(rpath)
            if not rd.empty:
                rd["i_evalue"] = pd.to_numeric(rd.i_evalue, errors="coerce")
                rcov = domain_coverage(rd)
                for _, fam in rescue_fams.iterrows():
                    sel = rcov[rcov.pfam_id.isin(fam.pfam_ids)]
                    for prot, g in sel.groupby("protein"):
                        if (fam.family, prot) in called:
                            continue
                        sub = rd[(rd.target == prot) & (rd.pfam_id.isin(fam.pfam_ids))]
                        cands.append({"fam": fam, "protein": prot,
                                      "prot_len": int(g.prot_len.iloc[0]),
                                      "hmm_cov": float(g.hmm_cov.max()),
                                      "best_evalue": float(sub.i_evalue.min()),
                                      "copies": _copies(sub)})

        # one pass over the proteome: size, short-protein capacity,
        # composition screen, and sequences for rescue candidates
        need = {x["protein"] for x in cands}
        seqs = {}
        faa = faa_of.get(sp)
        if faa:
            lens = []
            for h, seq in read_fasta(faa):
                pid = h.split()[0]
                L = len(seq.rstrip("*"))
                lens.append(L)
                if pid in need:
                    seqs[pid] = seq
                if has_mt and L <= 120 and _is_mt_like(seq):
                    f_ = mt_features(seq)
                    comp_rows.append((sp, pid, f_["length"], f_["n_cys"], f_["cys_fraction"],
                                      f_["n_motifs"], f_["n_aromatic"]))
            prot_totals[sp] = len(lens)
            lens.sort()
            stat_rows.append((sp, len(lens), sum(1 for x in lens if x <= 60),
                              sum(1 for x in lens if x <= 100),
                              lens[len(lens) // 2] if lens else 0,
                              lens[0] if lens else 0))

        for x in cands:
            fam = x["fam"]
            s_ = seqs.get(x["protein"])
            vals = {"len": x["prot_len"], "copies": x["copies"],
                    "cys": (s_.upper().count("C") / max(len(s_), 1)) if s_ else None}
            rules = _parse_validate(fam.validate, fam.family)
            ok, why = (x["hmm_cov"] >= fam.min_cov and x["prot_len"] >= fam.min_len,
                       "below map min_cov/min_len")
            if ok:
                ok, why = _check_rules(rules, vals)
            resc_rows.append((sp, fam.family, x["protein"], x["prot_len"],
                              round(x["hmm_cov"], 3), x["best_evalue"], x["copies"],
                              None if vals["cys"] is None else round(vals["cys"], 3),
                              "PASS" if ok else "FAIL", why))

        c.to_csv(os.path.join(calls_dir, f"{sp}.tsv"), sep="\t", index=False)
        all_calls.append(c)
        log(f"{sp}: {c.protein.nunique()} defensome proteins in {c.family.nunique()} families"
            + (f"; {sum(1 for r in resc_rows if r[0]==sp and r[8]=='PASS')} rescued"
               if len(rescue_fams) else ""))

    calls = pd.concat(all_calls, ignore_index=True)
    calls.to_csv(os.path.join(a.out, "gene_calls_all.tsv"), sep="\t", index=False)
    samples = [re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(p)) for p in files]
    counts = (calls.groupby(["species", "family"]).protein.nunique()
              .unstack(fill_value=0).reindex(index=samples, columns=dm.family, fill_value=0))
    counts.index.name = "species"
    counts.to_csv(os.path.join(a.out, "counts.tsv"), sep="\t")

    rc = pd.DataFrame(resc_rows, columns=["species", "family", "protein", "prot_len",
                                          "hmm_cov", "best_evalue", "copies",
                                          "cys_fraction", "validation", "reason"])
    rc.to_csv(os.path.join(a.out, "rescue_calls.tsv"), sep="\t", index=False)
    ok = rc[rc.validation == "PASS"]
    add = (ok.groupby(["species", "family"]).protein.nunique().unstack(fill_value=0)
             .reindex(index=counts.index, columns=counts.columns, fill_value=0))
    (counts + add).to_csv(os.path.join(a.out, "counts_rescued.tsv"), sep="\t")
    if len(rescue_fams):
        for f in rescue_fams.family:
            n_ga, n_r = int(counts[f].sum()), int(add[f].sum()) if f in add else 0
            log(f"{f}: {n_ga} at gathering threshold, {n_r} more rescued below it "
                f"({int((rc.family == f).sum())} candidates examined)")

    pd.DataFrame(comp_rows, columns=["species", "protein", "length", "n_cys", "cys_fraction",
                                     "n_cys_motifs", "n_aromatic"]) \
      .to_csv(os.path.join(a.out, "composition_candidates.tsv"), sep="\t", index=False)
    if stat_rows:
        pd.DataFrame(stat_rows, columns=["species", "n_proteins", "n_le60aa", "n_le100aa",
                                         "median_length", "min_length"]) \
          .to_csv(os.path.join(a.out, "proteome_stats.tsv"), sep="\t", index=False)

    if prot_totals:
        tot = pd.Series(prot_totals).reindex(counts.index)
        norm = counts.div(tot, axis=0) * 10000
        norm.round(3).to_csv(os.path.join(a.out, "counts_per10k.tsv"), sep="\t")
        pd.DataFrame({"species": tot.index, "n_proteins": tot.values}).to_csv(
            os.path.join(a.out, "proteome_sizes.tsv"), sep="\t", index=False)
        log("wrote counts.tsv, counts_rescued.tsv and counts_per10k.tsv")
    else:
        log("wrote counts.tsv and counts_rescued.tsv (no --proteomes, so no normalisation "
            "and rescue rules needing sequence cannot be checked)")
    log(f"annotate done: {counts.shape[0]} samples x {counts.shape[1]} families")


# --------------------------------------------------------------------------- #
def cmd_report(a):
    need_pandas()
    counts = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    dm = read_map(a.map)
    core = list(dm.loc[dm.tier == "CORE", "family"])

    summary = pd.DataFrame({
        "family": counts.columns,
        "tier": dm.set_index("family").tier.reindex(counts.columns).values,
        "min": counts.min(), "median": counts.median(),
        "max": counts.max(), "mean": counts.mean().round(2),
        "cv": (counts.std() / counts.mean().replace(0, np.nan)).round(3),
        "n_species_absent": (counts == 0).sum(),
    }).reset_index(drop=True).sort_values("cv", ascending=False)
    summary.to_csv(os.path.join(a.out, "family_summary.tsv"), sep="\t", index=False)
    print("\n=== Most variable families (CORE tier) ===")
    print(summary[summary.tier == "CORE"].head(10).to_string(index=False))

    if a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        if a.group_by not in md.columns:
            sys.exit(f"ERROR: column '{a.group_by}' not in {a.metadata}.\n"
                     f"Available columns: {list(md.columns)}")
        # Find whichever column actually matches the proteome names. Metadata
        # often stores "Genus species" while the file is "Genus_species.faa".
        want, key, best = set(counts.index), None, 0
        for c in md.columns:
            v = md[c].astype(str).str.strip().str.replace(" ", "_", regex=False)
            n = len(want & set(v))
            if n > best:
                key, best = c, n
        if key is None or best == 0:
            sys.exit(f"ERROR: no column in {a.metadata} matches the proteome names.\n"
                     f"Expected values like: {sorted(want)[:3]}\n"
                     f"Columns checked: {list(md.columns)}")
        log(f"metadata key column: '{key}' ({best}/{len(want)} species matched)")
        md = md.assign(**{key: md[key].astype(str).str.strip()
                          .str.replace(" ", "_", regex=False)}).set_index(key)
        use = norm if norm is not None else counts
        g = md[a.group_by].reindex(use.index)
        miss = g[g.isna()].index.tolist()
        if miss:
            log(f"WARNING: no metadata for {len(miss)} species: {miss[:5]}")
        gm = use.groupby(g).mean().round(2).T
        gm["n_groups"] = gm.notna().sum(axis=1)
        gm.to_csv(os.path.join(a.out, f"group_means_{a.group_by}.tsv"), sep="\t")
        print(f"\n=== Mean copies per 10k {size_unit(a.out)} by {a.group_by} (CORE) ===")
        print(gm.loc[[f for f in core if f in gm.index]].to_string())
        try:
            import warnings
            from scipy import stats
            rows = []
            for fam in use.columns:
                grp = [use.loc[g == lv, fam].values for lv in g.dropna().unique()]
                grp = [x for x in grp if len(x) > 1]
                if len(grp) >= 2 and any(len(set(x)) > 1 for x in grp):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        h, p = stats.kruskal(*grp)
                    if np.isfinite(p):
                        rows.append((fam, h, p))
            kw = pd.DataFrame(rows, columns=["family", "H", "p"]).sort_values("p")
            kw["p_bh"] = kw.p * len(kw) / kw.p.rank()   # Benjamini-Hochberg
            kw["p_bh"] = kw.p_bh[::-1].cummin()[::-1].clip(upper=1)
            kw.to_csv(os.path.join(a.out, f"kruskal_{a.group_by}.tsv"), sep="\t", index=False)
            if kw.empty:
                print(f"\nNo group test run: '{a.group_by}' needs at least two "
                      "groups with more than one species each.")
            else:
                print(f"\nKruskal-Wallis by {a.group_by}, top 5 (see kruskal_*.tsv):")
                print(kw.head(5).to_string(index=False))
            print("\nNOTE: species are not independent. These tests ignore phylogeny\n"
                  "      and are a screen, not a result.")
        except ImportError:
            log("scipy not available, skipping group tests")

    meta_df = None
    if a.metadata and a.group_by:
        meta_df = md
    make_figures(a.out, counts, norm, dm, meta_df, a.group_by)
    log("report done")


# =========================================================================== #
VERDICTS = {
    "ACCESSION_NOT_IN_DATABASE":
        "the HMM is not in the database that was searched; check the accession "
        "and the Pfam release in MANIFEST.tsv",
    "FILTERED_BY_MAP":
        "hits passed Pfam's threshold but every one failed this family's min_cov, "
        "min_len or ALL rule; inspect them before loosening the map",
    "FOUND_BELOW_GA":
        "no hits at Pfam's gathering threshold, but validated hits below it; "
        "reported in counts_rescued.tsv and rescue_calls.tsv, not in counts.tsv",
    "BELOW_GA_FAILED_VALIDATION":
        "below-threshold hits exist but none passed the map's validate rules; "
        "read rescue_calls.tsv to see why",
    "RECOVERED_FROM_TRANSCRIPTS":
        "absent from the proteomes, but MT-like ORFs were recovered from the transcripts by "
        "`short-orfs`; see short_orf_mt.tsv. The proteomes lost them to the ORF caller's length cut-off",
    "INPUT_LACKS_SHORT_PROTEINS":
        "a short family, and most proteomes contain almost no proteins of that "
        "length; the gene predictor's minimum length has probably removed them "
        "(TransDecoder's default is 100 aa), so absence cannot be concluded",
    "COMPOSITION_CANDIDATES_ONLY":
        "no HMM evidence, but proteins with metallothionein-like composition; "
        "see composition_candidates.tsv; heuristic, verify by alignment",
    "NOT_DETECTED":
        "searched with the right HMM, below threshold too, nothing found; "
        "consistent with genuine absence or extreme divergence",
    "NOT_DETECTED_AT_GA_ONLY":
        "nothing at Pfam's threshold, and the below-threshold rescue was not run "
        "for this family (rescue=no in the map, or scan --no-rescue); absence "
        "has not been tested as thoroughly as it could be",
}


def _manifest_for(out):
    p = os.path.join(out, "scan_db.txt")
    if not os.path.exists(p):
        return None
    db = open(p).read().strip()
    m = os.path.join(os.path.dirname(db), "MANIFEST.tsv")
    return pd.read_csv(m, sep="\t") if os.path.exists(m) else None


def diagnose_zero_families(a, c):
    """Turn every family that came back empty into a labelled verdict.

    A zero can mean the HMM was never searched, that a map rule filtered
    everything, that real copies score just below Pfam's threshold, that the
    input proteomes cannot contain proteins that short, or genuine absence.
    Those need different responses, and a count table cannot tell them apart."""
    dmap = read_map(a.map).set_index("family")
    dead = [f for f in c.columns if f in dmap.index and c[f].sum() == 0]
    if not dead:
        return
    man = _manifest_for(a.out)
    in_db = set(man[man.found == "yes"].pfam_id) if man is not None else None
    leng = dict(zip(man.pfam_id, man.length)) if man is not None else {}

    def raw_hits(sub, accs):
        n = {x: 0 for x in accs}
        for p in glob.glob(os.path.join(a.out, sub, "*.domtblout*")):
            with _open_text(p) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split(None, 6)
                    if len(f) > 4 and f[4].split(".")[0] in n:
                        n[f[4].split(".")[0]] += 1
        return n
    accs = {x for f in dead for x in dmap.loc[f, "pfam_ids"]}
    ga = raw_hits("hmmsearch", accs)
    rs = raw_hits("rescue", accs)
    rescue_ran = bool(glob.glob(os.path.join(a.out, "rescue", "*.domtblout*")))
    rc = _read_tsv(os.path.join(a.out, "rescue_calls.tsv"))
    cc = _read_tsv(os.path.join(a.out, "composition_candidates.tsv"))
    ps = _read_tsv(os.path.join(a.out, "proteome_stats.tsv"))
    so = _read_tsv(os.path.join(a.out, "short_orf_summary.tsv"))
    n_so_total = int(so.mt_like_genes.sum()) if so is not None and len(so) else 0

    rows = []
    for f in dead:
        ids = dmap.loc[f, "pfam_ids"]
        n_ga = sum(ga[x] for x in ids)
        n_rs = sum(rs[x] for x in ids)
        n_pass = int(((rc.family == f) & (rc.validation == "PASS")).sum()) if rc is not None and len(rc) else 0
        n_comp = len(cc) if (cc is not None and f == "Metallothionein") else 0
        L = min([int(leng[x]) for x in ids if x in leng and str(leng[x]).isdigit()] or [999])
        short = L < 60 or int(dmap.loc[f, "min_len"]) < 60
        lacking = None
        if ps is not None and len(ps):
            lacking = int((ps.n_le60aa <= 5).sum())
        missing = [x for x in ids if in_db is not None and x not in in_db]
        if missing:
            v = "ACCESSION_NOT_IN_DATABASE"
        elif n_ga > 0:
            v = "FILTERED_BY_MAP"
        elif n_pass > 0:
            v = "FOUND_BELOW_GA"
        elif f == "Metallothionein" and n_so_total > 0:
            v = "RECOVERED_FROM_TRANSCRIPTS"
        elif n_rs > 0:
            v = "BELOW_GA_FAILED_VALIDATION"
        elif short and lacking is not None and lacking >= max(1, len(ps) // 2):
            v = "INPUT_LACKS_SHORT_PROTEINS"
        elif n_comp > 0:
            v = "COMPOSITION_CANDIDATES_ONLY"
        elif not (rescue_ran and dmap.loc[f, "rescue"] == "yes"):
            v = "NOT_DETECTED_AT_GA_ONLY"
        else:
            v = "NOT_DETECTED"
        rows.append({"family": f, "pfam_ids": ",".join(ids), "verdict": v,
                     "in_database": "unknown" if in_db is None else ("no" if missing else "yes"),
                     "hmm_length": "" if L == 999 else L,
                     "raw_hits_at_GA": n_ga, "raw_hits_below_GA": n_rs,
                     "rescued_validated": n_pass, "composition_candidates": n_comp,
                     "short_orf_candidates": n_so_total if f == "Metallothionein" else 0,
                     "proteomes_with_<=5_proteins_under_60aa":
                         "" if lacking is None else f"{lacking}/{len(ps)}",
                     "what_it_means": VERDICTS[v]})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "qc_zero_families.tsv"), sep="\t", index=False)
    print("\n=== Families with no calls at Pfam's gathering threshold ===")
    for r in rows:
        print(f"  {r['family']:18s} {r['pfam_ids']:10s} {r['verdict']}")
        print(f"  {'':18s} {r['what_it_means']}")
    if man is None:
        print("  (no MANIFEST.tsv found next to the scanned database; run `setup` so "
              "missing accessions can be ruled out)")




def cmd_qc(a):
    """Annotation-quality report from low-copy control families.

    Two separate questions, reported separately:

    (1) Is a FAMILY behaving? Compares the median count across species against
        the literature expectation. If catalase has a median of 3 in every
        species, the map rule is catching something extra; that is a problem
        with the map, not with any one genome.

    (2) Is a SPECIES behaving? Compares each species against the observed
        median for that family, not against my prior. A genome with seven times
        the median catalase count has duplicated gene models. This is
        self-calibrating, so it stays valid on taxa where my priors do not.
    """
    need_pandas()
    check_fresh(a.out, read_map(a.map), a.map)
    c = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    diagnose_zero_families(a, c)
    sizes_f = os.path.join(a.out, "proteome_sizes.tsv")
    sizes = (pd.read_csv(sizes_f, sep="\t").set_index("species").n_proteins
             if os.path.exists(sizes_f) else None)
    # literature expectation for insect genomes, used only for question (1)
    # Literature expectation for insect genomes. Where the observed median
    # stays high after tightening the thresholds, the domain is simply shared
    # by a wider set of proteins than the family name implies, and the family
    # needs phylogeny to subdivide rather than a stricter cutoff.
    expect = {"GCL": 1, "Catalase": 1, "Hsp90": 2, "SOD_Fe": 1, "SOD_CuZn": 3,
              "NR_xeno": 21, "GST_MAPEG": 4, "BTB_BACK_Kelch": 7, "Keap1": 1}
    ctrl = [k for k in expect if k in c.columns]
    if not ctrl:
        sys.exit("ERROR: no low-copy control families found in counts.tsv")

    fam = pd.DataFrame({
        "expected": pd.Series({k: expect[k] for k in ctrl}),
        "observed_median": c[ctrl].median(),
        "observed_min": c[ctrl].min(), "observed_max": c[ctrl].max(),
    })
    fam["median_vs_expected"] = (fam.observed_median / fam.expected).round(1)
    # A median of zero has several possible causes (never searched, input too
    # short, genuinely absent) and the zero-family diagnosis above names the
    # right one. Calling every zero MAP_RULE_TOO_STRICT sent users to loosen a
    # rule for families whose HMM had never been in the database at all.
    zdiag = _read_tsv(os.path.join(a.out, "qc_zero_families.tsv"))
    zv = dict(zip(zdiag.family, zdiag.verdict)) if zdiag is not None and len(zdiag) else {}
    def _verdict(f, r):
        if r.observed_median == 0:
            return ("NOT_SEARCHED" if zv.get(f) == "ACCESSION_NOT_IN_DATABASE"
                    else "ZERO_SEE_DIAGNOSIS")
        if r.median_vs_expected >= 2:
            return "MAP_RULE_TOO_LOOSE"
        if r.median_vs_expected <= 0.5:
            return "MAP_RULE_TOO_STRICT"
        return "ok"
    fam["verdict"] = [_verdict(f, r) for f, r in fam.iterrows()]
    fam.sort_values("median_vs_expected", ascending=False) \
       .to_csv(os.path.join(a.out, "qc_families.tsv"), sep="\t")
    print("\n=== (1) Are the control FAMILIES behaving? ===")
    print(fam.sort_values("median_vs_expected", ascending=False).to_string())
    rule = fam[fam.verdict.isin(["MAP_RULE_TOO_LOOSE", "MAP_RULE_TOO_STRICT"])]
    if len(rule):
        print(f"\n-> {list(rule.index)}: counts are present but off the expectation, so the")
        print("   map rule is the likely cause. Adjust min_cov/min_len or require a partner domain.")
    unsearched = fam[fam.verdict == "NOT_SEARCHED"]
    if len(unsearched):
        print(f"\n-> {list(unsearched.index)}: never searched; the HMM is not in the database.")
        print("   Fix the database (run `setup` against a full Pfam-A), not the map.")
    zero = fam[fam.verdict == "ZERO_SEE_DIAGNOSIS"]
    if len(zero):
        print(f"\n-> {list(zero.index)}: zero everywhere; the zero-family table above says why.")

    # question (2): species outliers against the observed median
    base = c[ctrl].median().replace(0, np.nan)
    usable = [k for k in ctrl if pd.notna(base[k])]
    if not usable:
        print("\n=== (2) Species outliers: not assessable ===")
        print("Every control family has a median of zero across samples, so no")
        print("species can be compared against it. Either the controls are genuinely")
        print("absent from these proteomes or the map rules for them are too strict;")
        print("the table under (1) shows which.")
        log("qc done")
        return
    ratio = c[usable].div(base[usable], axis=1)
    # pandas >= 3 raises on idxmax over an all-NA row; compute it only where
    # at least one control family is informative for that sample
    has = ratio.notna().any(axis=1)
    worst = pd.Series("", index=ratio.index, dtype=object)
    if has.any():
        worst[has] = ratio[has].idxmax(axis=1)
    qc = pd.DataFrame({
        "n_proteins": sizes.reindex(c.index) if sizes is not None else np.nan,
        "max_ratio": ratio.max(axis=1).round(1),
        "worst_family": worst,
        "n_families_inflated": (ratio >= 2).sum(axis=1),
        "n_families_missing": (c[ctrl] == 0).sum(axis=1),
    })
    for k in ctrl:
        qc[f"n_{k}"] = c[k]
    qc["flag"] = np.where((qc.max_ratio >= 3) | (qc.n_families_inflated >= 2), "INFLATED",
                  np.where(qc.n_families_missing >= 2, "SPARSE", "ok"))
    qc = qc.sort_values("max_ratio", ascending=False)
    qc.to_csv(os.path.join(a.out, "qc_annotation.tsv"), sep="\t")
    print("\n=== (2) Are any SPECIES outliers? (ratio to observed median) ===")
    print(qc.head(10).to_string())
    n_bad = (qc.flag != "ok").sum()
    print(f"\n{n_bad} of {len(qc)} species flagged. Use max_ratio as an "
          "annotation-quality covariate,\nor exclude the flagged species. "
          "This replaces a BUSCO duplication run for most purposes.")
    log("qc done")


# =========================================================================== #
# =========================================================================== #
#  Figure style, shared by every figure the tool writes
# =========================================================================== #
# Same colours as the dashboard, so a method means one colour everywhere.
METHOD_COLOURS = {"DToL": "#1d3557", "genome_guided": "#457b9d", "denovo": "#e63946"}
_METHOD_RANK = {"DToL": 0, "genome": 0, "genome_guided": 1, "denovo": 2}


def method_order(methods):
    """Genome annotation first, then genome-guided, then de novo: the order of
    decreasing reliance on a genome, so lines between them read as a gradient
    rather than zig-zagging alphabetically."""
    return sorted(methods, key=lambda m: (_METHOD_RANK.get(m, 9), str(m)))


def method_colour(m, i=0):
    fallback = ["#6a994e", "#bc6c25", "#7b2cbf", "#8d99ae"]
    return METHOD_COLOURS.get(m, fallback[i % len(fallback)])


def _pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": 9,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e9ecef", "grid.linewidth": .6,
        "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        # editable text in Illustrator/Inkscape rather than outlined glyphs
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": "tight", "savefig.dpi": 300,
    })


class FigureWriter:
    """Every figure as PDF (vector, for the thesis) and PNG (300 dpi, for
    slides), plus one combined PDF and a captions file so each figure says
    what it shows without the reader having to open this code."""
    def __init__(self, fig_dir, summary_name):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        self.plt, self.dir = plt, fig_dir
        os.makedirs(fig_dir, exist_ok=True)
        self.pages = PdfPages(os.path.join(fig_dir, summary_name))
        self.captions = []

    def save(self, fig, name, caption=""):
        base = os.path.join(self.dir, name)
        fig.savefig(base + ".pdf")
        fig.savefig(base + ".png")
        self.pages.savefig(fig)
        self.plt.close(fig)
        if caption:
            self.captions.append((name, caption))

    def close(self, caption_file="FIGURES.md"):
        self.pages.close()
        with open(os.path.join(self.dir, caption_file), "a") as fh:
            for n, c in self.captions:
                fh.write(f"**{n}**  {c}\n\n")



def _split_families(dm, mat):
    """CORE families that can be z-scored, and those that cannot.

    An all-zero or constant family has no z-score. Dropping it silently is
    wrong: 'GCL is single copy in every species' and 'MATE was never detected'
    are both findings and both belong in the output.
    """
    core = [f for f in dm.loc[dm.tier == "CORE", "family"] if f in mat.columns]
    var, absent, const = [], [], []
    for f in core:
        if mat[f].sum() == 0:
            absent.append(f)
        elif mat[f].std() == 0:
            const.append(f)
        else:
            var.append(f)
    return core, var, absent, const


def make_figures(out, counts, norm, dm, meta=None, group=None):
    np.seterr(divide="ignore", invalid="ignore")
    try:
        _pub_style()
    except Exception:
        pass
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        log("matplotlib not available, skipping figures (pip install --user matplotlib)")
        return
    fig_dir = os.path.join(out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    use = norm if norm is not None else counts
    unit = size_unit(out)
    lab = f"copies per 10k {unit}" if norm is not None else "copies"
    cats = dm.set_index("family").category

    # every count figure is produced twice: normalised and raw
    variants = [("", use, lab)]
    if norm is not None:
        variants.append(("_raw", counts, "raw copies"))

    core, var, _, _ = _split_families(dm, use)
    # "absent" and "constant" are statements about the genomes, so judge them on
    # RAW counts. GCL is exactly 1 in every species but varies after dividing by
    # proteome size, and that single-copy result is the calibration finding.
    _, _, absent, const = _split_families(dm, counts)
    if absent or const:
        log(f"invariant families held out of z-scored panels: "
            f"absent={absent} constant={const} (see 00_invariant_families.pdf)")

    lv, cmap = [], {}
    if meta is not None and group is not None and group in meta.columns:
        lv = sorted(meta[group].dropna().unique())
        cmap = dict(zip(lv, plt.cm.viridis(np.linspace(0, .85, len(lv)))))

    if os.path.exists(os.path.join(fig_dir, "FIGURES.md")):
        os.remove(os.path.join(fig_dir, "FIGURES.md"))
    W = FigureWriter(fig_dir, "report_all_figures.pdf")
    _caps = {'00_invariant': 'Families that cannot be z-scored because they are absent or constant in every sample, with their raw counts; judged on raw counts.', '01_heatmap': 'Samples by CORE family, colour is the z-score of each family across samples. Colour is relative to the other samples in this run, not an absolute amount.', '02_total': 'Total CORE defensome per sample.', '03_pca': 'PCA of z-scored CORE families: sample scores (left) and family loadings (right). Read the loadings to see which families drive each axis.', '04_category': 'Share of the CORE defensome in each functional tier, per sample.', '05_family': 'Pearson correlation between families across samples: blocks of red are families that expand together.', '06_normalisation': 'Evidence that dividing by proteome size removed the size confound: correlation with proteome size before (raw) and after (per 10k).', '07_family': 'Every CORE family by group, boxes with points, n per group in the axis labels.', '08_domain': 'Share of called proteins whose required domains are complete, truncated (PARTIAL) or unusually short (FRAGMENT), per family.'}

    def save(fig, name):
        fig.tight_layout()
        stem = name[:-4] if name.endswith(".pdf") else name
        cap = next((v for k, v in _caps.items() if stem.startswith(k)), "")
        if "_raw" in stem:
            cap = (cap + " Raw counts.").strip()
        W.save(fig, stem, cap)

    # ---- 00. the families a heatmap cannot show ----------------------------
    if absent or const:
        show = absent + const
        fig, axes = plt.subplots(1, len(show), figsize=(3.6*len(show), 4.6),
                                 squeeze=False)
        for ax, f in zip(axes[0], show):
            vals = counts[f]
            ax.bar(range(len(vals)), vals.values,
                   color="#adb5bd" if f in absent else "#2a6f97")
            ax.set_title(f, fontsize=11)
            ax.set_xticks([]); ax.set_xlabel(f"{len(vals)} species", fontsize=8)
            ax.set_ylabel("raw copies", fontsize=9)
            if f in absent:
                ax.set_ylim(0, 1)
                ax.text(.5, .55, "never detected\nin any species", ha="center",
                        va="center", transform=ax.transAxes, fontsize=9, color="#c1121f")
                ax.text(.5, .3, "run `qc` to see whether the\nHMM produced raw hits at all",
                        ha="center", va="center", transform=ax.transAxes, fontsize=7.5,
                        color="#666")
            else:
                ax.set_ylim(0, max(2, vals.max()*1.6))
                ax.text(.5, .8, f"constant at {int(vals.iloc[0])}\nin all species",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color="#2a6f97")
                ax.text(.5, .6, "single-copy control:\nthe pipeline is calibrated",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=7.5, color="#666")
        fig.suptitle("Families excluded from the z-scored panels, and why", y=1.02)
        save(fig, "00_invariant_families.pdf")

    # ---- 01. heatmaps, normalised and raw ----------------------------------
    for suffix, mat_, lab_ in variants:
        _, v_, _, _ = _split_families(dm, mat_)
        mm = mat_[v_]
        zz = ((mm - mm.mean()) / mm.std().replace(0, np.nan)).fillna(0)
        os_ = list(zz.index[np.argsort(zz.values @ zz.values.sum(0))])
        of_ = list(zz.columns[np.argsort(zz.values.sum(0))])
        fig, ax = plt.subplots(figsize=(max(9, .40*len(v_)), max(7, .30*len(mm))))
        im = ax.imshow(zz.loc[os_, of_].values, aspect="auto", cmap="RdBu_r",
                       vmin=-3, vmax=3)
        ax.set_xticks(range(len(of_))); ax.set_xticklabels(of_, rotation=90, fontsize=8)
        ax.set_yticks(range(len(os_)))
        ax.set_yticklabels([s.replace("_", " ") for s in os_], fontsize=7, style="italic")
        if lv:
            for t_, s_ in zip(ax.get_yticklabels(), os_):
                t_.set_color(cmap.get(meta[group].get(s_), "#222"))
        fig.colorbar(im, ax=ax, label="z-score", shrink=.6)
        note = ""
        if absent or const:
            note = "\nnot shown: " + ", ".join(absent + const) + " (see 00_invariant_families.pdf)"
        ax.set_title(f"Defensome, CORE families ({lab_}, z-scored per family){note}",
                     fontsize=11)
        save(fig, f"01_heatmap_core{suffix}.pdf")

    # ---- 02. totals, normalised and raw ------------------------------------
    for suffix, mat_, lab_ in variants:
        t_ = mat_[var].sum(axis=1).sort_values()
        cc = [cmap.get(meta[group].get(s_), "#999999") for s_ in t_.index] if lv else "#2a6f97"
        fig, ax = plt.subplots(figsize=(7.5, max(6, .26*len(t_))))
        ax.barh(range(len(t_)), t_.values, color=cc)
        ax.set_yticks(range(len(t_)))
        ax.set_yticklabels([s.replace("_", " ") for s in t_.index], fontsize=7.5,
                           style="italic")
        ax.set_xlabel(f"total CORE defensome, {lab_}")
        if lv:
            ax.legend(handles=[Patch(color=cmap[k], label=str(k)) for k in lv],
                      fontsize=8, title=group)
        save(fig, f"02_total_defensome{suffix}.pdf")

    # ---- 03. PCA -----------------------------------------------------------
    # needs at least two varying families and three samples, or there is no
    # second axis to plot; small or narrow runs skip it rather than crash
    if len(var) >= 2 and len(use) >= 3:
        m = use[var]
        z = ((m - m.mean()) / m.std().replace(0, np.nan)).fillna(0)
        x = z.values - z.values.mean(0)
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        pc = u[:, :2] * s[:2]; ev = (s**2 / (s**2).sum() * 100)[:2]
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6.5))
        if lv:
            for k in lv:
                i = [j for j, sp in enumerate(z.index) if meta[group].get(sp) == k]
                ax.scatter(pc[i, 0], pc[i, 1], s=42, color=cmap[k], label=str(k),
                           edgecolor="white", linewidth=.6)
            ax.legend(fontsize=9, title=group)
        else:
            ax.scatter(pc[:, 0], pc[:, 1], s=34, color="#2a6f97")
        for i, sp in enumerate(z.index):
            ax.annotate(sp.replace("_", " "), (pc[i, 0], pc[i, 1]), fontsize=5.5,
                        xytext=(4, 2), textcoords="offset points", style="italic")
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}% of variance)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.axhline(0, lw=.4, c="grey"); ax.axvline(0, lw=.4, c="grey")
        ax.set_title(f"Species scores ({lab}, {len(var)} CORE families z-scored)")
        ld = vt[:2].T * s[:2] / np.sqrt(len(z))
        ax2.axhline(0, lw=.5, c="grey"); ax2.axvline(0, lw=.5, c="grey")
        for i, f in enumerate(z.columns):
            ax2.annotate("", xy=(ld[i, 0], ld[i, 1]), xytext=(0, 0),
                         arrowprops=dict(arrowstyle="->", color="#c1121f", lw=.7))
            ax2.annotate(f, (ld[i, 0], ld[i, 1]), fontsize=7)
        ax2.set_xlabel("PC1 loading"); ax2.set_ylabel("PC2 loading")
        ax2.set_title("Family loadings: which families drive the separation")
        save(fig, "03_pca.pdf")

    else:
        log(f"PCA skipped: needs >=2 varying CORE families and >=3 samples "
            f"(have {len(var)} and {len(use)})")

    # ---- 04. category composition, normalised and raw ----------------------
    for suffix, mat_, lab_ in variants:
        cs = mat_[var].T.groupby(cats.reindex(var).values).sum().T
        cs = cs.loc[cs.sum(axis=1).sort_values().index]
        frac = cs.div(cs.sum(axis=1), axis=0)
        fig, ax = plt.subplots(figsize=(9.5, max(6, .26*len(frac))))
        left = np.zeros(len(frac))
        for c_ in frac.columns:
            ax.barh(range(len(frac)), frac[c_].values, left=left, label=c_)
            left += frac[c_].values
        ax.set_yticks(range(len(frac)))
        ax.set_yticklabels([s.replace("_", " ") for s in frac.index], fontsize=7.5,
                           style="italic")
        ax.set_xlabel(f"fraction of CORE defensome ({lab_})")
        ax.legend(fontsize=8, ncol=3, loc="lower right")
        save(fig, f"04_category_composition{suffix}.pdf")

    # ---- 05. family co-variation ------------------------------------------
    if len(var) >= 2 and len(use) >= 3:
        with np.errstate(divide="ignore", invalid="ignore"):
            cm = np.nan_to_num(np.corrcoef(z.values.T))
        o = np.argsort(cm.sum(0))
        fig, ax = plt.subplots(figsize=(max(8, .36*len(var)), max(7, .36*len(var))))
        im = ax.imshow(cm[np.ix_(o, o)], cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(var))); ax.set_xticklabels([var[i] for i in o], rotation=90, fontsize=8)
        ax.set_yticks(range(len(var))); ax.set_yticklabels([var[i] for i in o], fontsize=8)
        fig.colorbar(im, ax=ax, label="Pearson r", shrink=.7)
        ax.set_title(f"Do defensome families expand together? ({lab})")
        save(fig, "05_family_correlation.pdf")


    # ---- 06. normalisation QC ----------------------------------------------
    sizes_f = os.path.join(out, "proteome_sizes.tsv")
    if os.path.exists(sizes_f) and norm is not None:
        sz = pd.read_csv(sizes_f, sep="\t").set_index("species").n_proteins.reindex(counts.index)
        def _r(xx, yy):
            return 0.0 if np.std(xx) == 0 or np.std(yy) == 0 else float(np.corrcoef(xx, yy)[0, 1])
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
        axes[0].scatter(sz, counts[var].sum(axis=1), s=30, c="#c1121f")
        axes[0].set_xlabel(f"{unit} in proteome"); axes[0].set_ylabel("raw CORE defensome")
        axes[0].set_title(f"raw counts track proteome size (r={_r(counts[var].sum(axis=1), sz):.2f})")
        axes[1].scatter(sz, norm[var].sum(axis=1), s=30, c="#2a6f97")
        axes[1].set_xlabel(f"{unit} in proteome"); axes[1].set_ylabel(f"CORE per 10k {unit}")
        axes[1].set_title(f"normalised counts should not (r={_r(norm[var].sum(axis=1), sz):.2f})")
        rr = sorted([(f, _r(counts[f], sz), _r(norm[f], sz)) for f in var], key=lambda t: -t[1])
        axes[2].barh(range(len(rr)), [t[1] for t in rr], .4, label="raw", color="#c1121f")
        axes[2].barh([i+.4 for i in range(len(rr))], [t[2] for t in rr], .4,
                     label="per 10k", color="#2a6f97")
        axes[2].set_yticks([i+.2 for i in range(len(rr))])
        axes[2].set_yticklabels([t[0] for t in rr], fontsize=7)
        axes[2].axvline(0, c="k", lw=.6); axes[2].set_xlabel(f"r with number of {unit}")
        axes[2].legend(fontsize=8)
        save(fig, "06_normalisation_qc.pdf")

    # ---- 07. every family by group, paginated, full labels -----------------
    if lv:
        g = meta[group].reindex(use.index)
        counts_per_group = {k: int((g == k).sum()) for k in lv}
        for suffix, mat_, lab_ in variants:
            sel = core                      # ALL core families, including invariant
            per_page, nc = 20, 4
            pages = [sel[i:i+per_page] for i in range(0, len(sel), per_page)]
            for pi, page in enumerate(pages, 1):
                nr = int(np.ceil(len(page)/nc))
                fig, axes = plt.subplots(nr, nc, figsize=(4.3*nc, 3.3*nr), squeeze=False)
                for ax_, f in zip(np.ravel(axes), page):
                    data = [mat_.loc[g == k, f].values for k in lv]
                    bp = ax_.boxplot(data, positions=range(1, len(lv)+1),
                                     widths=.55, showfliers=False, patch_artist=True)
                    for patch, k in zip(bp["boxes"], lv):
                        patch.set_facecolor(cmap[k]); patch.set_alpha(.35)
                    for i, dv in enumerate(data, 1):
                        ax_.scatter(np.random.normal(i, .07, len(dv)), dv, s=11,
                                    color=cmap[lv[i-1]], alpha=.85, zorder=3,
                                    edgecolor="white", linewidth=.3)
                    ax_.set_xticks(range(1, len(lv)+1))
                    ax_.set_xticklabels([f"{k}\n(n={counts_per_group[k]})" for k in lv],
                                        fontsize=7.5)
                    ax_.set_title(f, fontsize=10)
                    ax_.tick_params(axis="y", labelsize=7.5)
                    if mat_[f].sum() == 0:
                        ax_.text(.5, .5, "never detected", ha="center", va="center",
                                 transform=ax_.transAxes, fontsize=9, color="#c1121f")
                    elif mat_[f].std() == 0:
                        ax_.text(.5, .85, "constant", ha="center", va="center",
                                 transform=ax_.transAxes, fontsize=8, color="#2a6f97")
                for ax_ in np.ravel(axes)[len(page):]:
                    ax_.axis("off")
                fig.suptitle(f"CORE families by {group} ({lab_})  "
                             f"page {pi} of {len(pages)}", y=1.005, fontsize=13)
                save(fig, f"07_family_by_group{suffix}_p{pi}.pdf")

    # ---- 08. domain completeness -------------------------------------------
    comp_f = os.path.join(out, "domains", "completeness.tsv")
    if os.path.exists(comp_f):
        comp = pd.read_csv(comp_f, sep="\t")
        t = comp.groupby(["family", "status"]).size().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in t:
                t[s_] = 0
        n_calls = t.sum(axis=1)
        t = t[["COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"]].div(n_calls, axis=0)
        t = t.loc[t.COMPLETE.sort_values().index]
        fig, ax = plt.subplots(figsize=(10, max(6, .28*len(t))))
        left = np.zeros(len(t))
        for c_, col in zip(t.columns, ["#2a6f97", "#f4a261", "#e76f51", "#adb5bd"]):
            ax.barh(range(len(t)), t[c_].values, left=left, label=c_, color=col)
            left += t[c_].values
        ax.set_yticks(range(len(t)))
        ax.set_yticklabels([f"{f}  (n={int(n_calls[f])})" for f in t.index], fontsize=8)
        ax.set_xlabel("fraction of called proteins"); ax.set_xlim(0, 1)
        ax.legend(fontsize=8, ncol=4, loc="lower left")
        ax.set_title("Domain architecture completeness per family")
        save(fig, "08_domain_completeness.pdf")


    # ---- R1. what the below-threshold rescue added -----------------------
    rc_f = os.path.join(out, "rescue_calls.tsv")
    cr_f = os.path.join(out, "counts_rescued.tsv")
    if os.path.exists(rc_f) and os.path.exists(cr_f):
        rc = pd.read_csv(rc_f, sep="\t")
        if len(rc):
            cr = pd.read_csv(cr_f, sep="\t", index_col=0)
            fams_r = sorted(rc.family.unique())
            fig, axes = plt.subplots(1, len(fams_r), figsize=(4.4 * len(fams_r), max(3.2, .22 * len(counts) + 1.4)),
                                     squeeze=False)
            for ax, f in zip(axes[0], fams_r):
                ga = counts[f] if f in counts else pd.Series(0, index=counts.index)
                tot_r = cr[f].reindex(counts.index).fillna(0) if f in cr else ga
                order = tot_r.sort_values().index
                y = np.arange(len(order))
                ax.barh(y, ga.reindex(order).values, color="#1d3557", label="at gathering threshold")
                ax.barh(y, (tot_r - ga).reindex(order).values, left=ga.reindex(order).values,
                        color="#f4a261", label="rescued below it")
                ax.set_yticks(y); ax.set_yticklabels([s_.replace("_", " ") for s_ in order], fontsize=7)
                sub = rc[rc.family == f]
                ax.set_title(f"{f}\n{int((sub.validation=='PASS').sum())} rescued, "
                             f"{int((sub.validation=='FAIL').sum())} rejected", loc="left", fontsize=10)
                ax.set_xlabel("genes"); ax.grid(axis="y", visible=False)
            from matplotlib.patches import Patch
            fig.legend(handles=[Patch(color="#1d3557", label="at Pfam's gathering threshold"),
                                Patch(color="#f4a261", label="rescued below it, validated")],
                       loc="upper center", ncol=2, bbox_to_anchor=(.5, 1.06), fontsize=8.5)
            save(fig, "R1_rescued_families.pdf")
            W.captions[-1] = ("R1_rescued_families",
                "Genes called at Pfam's gathering threshold (dark) and the additional genes found "
                "by the below-threshold rescue that passed the map's validation rules (orange). "
                "Rescued genes are reported in counts_rescued.tsv only, never in counts.tsv.")

    # ---- D1. why each empty family is empty --------------------------------
    z_f = os.path.join(out, "qc_zero_families.tsv")
    if os.path.exists(z_f):
        z0 = pd.read_csv(z_f, sep="\t")
        if len(z0):
            vcol = {"ACCESSION_NOT_IN_DATABASE": "#6a040f", "FILTERED_BY_MAP": "#e9a13b",
                    "FOUND_BELOW_GA": "#2a9d8f", "BELOW_GA_FAILED_VALIDATION": "#f4a261",
                    "INPUT_LACKS_SHORT_PROTEINS": "#7b2cbf", "COMPOSITION_CANDIDATES_ONLY": "#8d99ae",
                    "NOT_DETECTED": "#adb5bd", "NOT_DETECTED_AT_GA_ONLY": "#ced4da"}
            z0 = z0.sort_values(["verdict", "family"])
            fig, ax = plt.subplots(figsize=(9.5, .38 * len(z0) + 1.3))
            for yi, r in enumerate(z0.itertuples(index=False)):
                ax.add_patch(plt.Rectangle((0, yi - .4), 1, .8, color=vcol.get(r.verdict, "#adb5bd")))
                ax.text(-.03, yi, f"{r.family}  ({r.pfam_ids})", ha="right", va="center", fontsize=8.5)
                label = (r.verdict.replace("_", " ").lower()
                         .replace(" ga", " GA").replace("map", "map rules"))
                light = r.verdict in ("NOT_DETECTED", "NOT_DETECTED_AT_GA_ONLY",
                                      "COMPOSITION_CANDIDATES_ONLY")
                # text is always dark enough to read; the swatch carries the colour
                ax.text(1.05, yi - .1, label, va="center", fontsize=9, fontweight="bold",
                        color="#495057" if light else vcol.get(r.verdict, "#495057"))
                detail = (f"GA hits {r.raw_hits_at_GA} · below-GA hits {r.raw_hits_below_GA} · "
                          f"rescued {r.rescued_validated}")
                ax.text(1.05, yi + .24, detail, va="center", fontsize=7.2, color="#6c757d")
            ax.set_xlim(-2.2, 5.2); ax.set_ylim(len(z0) - .5, -.7)
            ax.axis("off")
            ax.set_title("Why each empty family is empty", loc="left")
            save(fig, "D1_zero_family_verdicts.pdf")
            W.captions[-1] = ("D1_zero_family_verdicts",
                "Every family with no calls at Pfam's gathering threshold, with the pipeline's verdict "
                "on why: missing from the database, filtered by the map, found only below the "
                "threshold, absent because the input lacks proteins that short, or not detected. "
                "Full detail and advice in qc_zero_families.tsv.")

    W.close()
    log(f"figures written to {fig_dir} (PDF, PNG, report_all_figures.pdf, FIGURES.md)")


# =========================================================================== #
#  Sequence and domain utilities
# =========================================================================== #
def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def tool(name):
    """Absolute path to an external tool, or None.
    DEFENSOME_<NAME> (e.g. DEFENSOME_HMMSEARCH=/opt/hmmer/bin/hmmsearch) wins
    over PATH, for clusters where HMMER and Python cannot share a module set."""
    env = os.environ.get("DEFENSOME_" + name.upper().replace("-", "_"))
    if env and os.path.exists(env) and os.access(env, os.X_OK):
        return env
    return shutil.which(name)


def run(cmd, **kw):
    cmd = [str(c) for c in cmd]
    resolved = tool(cmd[0])
    if resolved:
        cmd[0] = resolved
    log("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


def all_domains(path):
    """Every Pfam domain on every protein, with coordinates. Unlike
    domain_coverage() this keeps each hit separately so architectures and
    per-domain completeness can be reconstructed."""
    op = gzip.open if path.endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 22)
            if len(f) < 22:
                continue
            rows.append((f[0], int(f[2]), f[3], f[4].split(".")[0], int(f[5]),
                         int(f[15]), int(f[16]), int(f[17]), int(f[18]), float(f[13])))
    return pd.DataFrame(rows, columns=["protein","prot_len","pfam_name","pfam_id",
                                       "hmm_len","hmm_from","hmm_to","ali_from",
                                       "ali_to","score"])


def check_fresh(out, dm, mapfile):
    """Stop if gene_calls_all.tsv was built with a different version of the map.

    Editing the map and rerunning a downstream command without redoing
    `annotate` silently mixes two family vocabularies. `annotate` is cheap
    because it reads the cached domtblout, so the fix is always the same.
    """
    calls_f = os.path.join(out, "gene_calls_all.tsv")
    if not os.path.exists(calls_f):
        sys.exit(f"ERROR: {calls_f} not found. Run `annotate` first.")
    fams = set(pd.read_csv(calls_f, sep="\t", usecols=["family"]).family.unique())
    orphan = sorted(fams - set(dm.family))
    if orphan:
        sys.exit(
            f"\nERROR: gene_calls_all.tsv is stale.\n"
            f"  Families present in the calls but not in {os.path.basename(mapfile)}: "
            f"{orphan}\n"
            f"  The map was edited after `annotate` last ran.\n\n"
            f"  Fix (seconds, reuses the cached hmmsearch output):\n"
            f"    python3 defensome.py annotate --proteomes <DIR> --out {out}\n")
    # mapfile may not exist on disk at all when the embedded default is in use
    if os.path.exists(mapfile) and os.path.getmtime(mapfile) > os.path.getmtime(calls_f):
        log(f"WARNING: {os.path.basename(mapfile)} is newer than gene_calls_all.tsv. "
            "Family names still match, but thresholds may have changed. "
            "Rerun `annotate` to be sure.")


def scan_files(out):
    fs = sorted(glob.glob(os.path.join(out, "hmmsearch", "*.domtblout*")))
    if not fs:
        sys.exit(f"ERROR: no domtblout files in {out}/hmmsearch. Run `scan` first.")
    return [(re.sub(r"\.domtblout(\.gz)?$", "", os.path.basename(f)), f) for f in fs]


# =========================================================================== #
def cmd_extract(a):
    """Per-family protein FASTA. Needed by `trees` and `cyp`."""
    need_pandas()
    check_fresh(a.out, read_map(a.map), a.map)
    calls_f = os.path.join(a.out, "gene_calls_all.tsv")
    calls = pd.read_csv(calls_f, sep="\t")
    if a.families:
        want = set(a.families.split(","))
        calls = calls[calls.family.isin(want)]
        unknown = want - set(calls.family)
        if unknown:
            log(f"WARNING: no calls for {sorted(unknown)}")
    fdir = os.path.join(a.out, "fasta"); os.makedirs(fdir, exist_ok=True)

    wanted = calls.groupby("species").protein.apply(set).to_dict()
    seqs = {}
    for sp, faa in species_list(a.proteomes):
        if sp not in wanted:
            continue
        keep = wanted[sp]
        for hdr, seq in read_fasta(faa):
            pid = hdr.split()[0]
            if pid in keep:
                seqs[(sp, pid)] = seq
    missing = sum(1 for sp, g in calls.groupby("species")
                  for p in g.protein.unique() if (sp, p) not in seqs)
    if missing:
        log(f"WARNING: {missing} called proteins not found in the proteome FASTAs "
            "(header/ID mismatch?)")

    n = 0
    for fam, g in calls.groupby("family"):
        path = os.path.join(fdir, f"{fam}.faa")
        with open(path, "w") as fh:
            for sp, prot in g[["species","protein"]].drop_duplicates().itertuples(index=False):
                s = seqs.get((sp, prot))
                if s:
                    fh.write(f">{sp}|{prot}\n{s}\n")
        n += 1
    log(f"extract done: {n} family FASTA files in {fdir}")


# =========================================================================== #
def cmd_domains(a):
    """Per-protein domain architecture and completeness.

    For every defensome protein this records each Pfam domain present, where it
    sits, and how much of the HMM it covers. A protein is COMPLETE when every
    domain its family requires is present above the family coverage threshold.
    Anything else is PARTIAL (domain present but truncated) or MISSING_DOMAIN.
    FRAGMENT flags proteins far shorter than the family norm, which is usually a
    split gene model rather than a real short paralogue."""
    need_pandas()
    dm_raw = read_map(a.map)
    check_fresh(a.out, dm_raw, a.map)
    dm = dm_raw.set_index("family")
    calls = pd.read_csv(os.path.join(a.out, "gene_calls_all.tsv"), sep="\t")
    ddir = os.path.join(a.out, "domains"); os.makedirs(ddir, exist_ok=True)

    want = calls.groupby("species").protein.apply(set).to_dict()
    dom_rows = []
    for sp, path in scan_files(a.out):
        if sp not in want:
            continue
        d = all_domains(path)
        d = d[d.protein.isin(want[sp])].copy()
        d.insert(0, "species", sp)
        dom_rows.append(d)
    dom = pd.concat(dom_rows, ignore_index=True)
    dom["dom_cov"] = ((dom.hmm_to - dom.hmm_from + 1) / dom.hmm_len).round(3)
    dom.to_csv(os.path.join(ddir, "domain_table.tsv"), sep="\t", index=False)
    log(f"domain_table.tsv: {len(dom)} domain hits on {dom.protein.nunique()} proteins")

    # Ordered architecture string per protein, e.g. GST_N-GST_C
    arch = (dom.sort_values(["species","protein","ali_from"])
              .groupby(["species","protein"])
              .agg(architecture=("pfam_name", lambda s: "-".join(s)),
                   n_domains=("pfam_name", "size"),
                   prot_len=("prot_len", "first")).reset_index())
    arch.to_csv(os.path.join(ddir, "architectures.tsv"), sep="\t", index=False)

    # Per-protein completeness against the family's required domains
    # Merge non-overlapping HMM segments per (protein, pfam) exactly as
    # annotate does. Many Pfam models match in two pieces (FMO's paired
    # Rossmann folds, ABC_membrane, thioredoxin); scoring on the best single
    # segment understates coverage and mislabels intact proteins as PARTIAL.
    best = {}
    for (sp_, pr_, pf_), g_ in dom.groupby(["species","protein","pfam_id"], sort=False):
        iv = sorted(zip(g_.hmm_from.astype(int), g_.hmm_to.astype(int)))
        cs = ce = None; tot = 0
        for s_, e_ in iv:
            if cs is None: cs, ce = s_, e_
            elif s_ <= ce + 1: ce = max(ce, e_)
            else: tot += ce - cs + 1; cs, ce = s_, e_
        if cs is not None: tot += ce - cs + 1
        hl = float(g_.hmm_len.iloc[0])
        best[(sp_, pr_, pf_)] = tot / hl if hl else 0.0
    med_len = calls.groupby("family").prot_len.median().to_dict()

    rows = []
    skipped = set()
    for r in calls.itertuples(index=False):
        if r.family not in dm.index:
            skipped.add(r.family); continue
        fam = dm.loc[r.family]
        req = list(fam.pfam_ids) if fam.rule == "ALL" else [
            p for p in fam.pfam_ids if (r.species, r.protein, p) in best]
        covs = {p: best.get((r.species, r.protein, p), 0.0) for p in req}
        n_ok = sum(1 for v in covs.values() if v >= fam.min_cov)
        if not req:
            status = "MISSING_DOMAIN"
        elif n_ok == len(req):
            status = "COMPLETE"
        elif any(v > 0 for v in covs.values()):
            status = "PARTIAL"
        else:
            status = "MISSING_DOMAIN"
        lr = r.prot_len / med_len.get(r.family, r.prot_len)
        if status == "COMPLETE" and lr < 0.5:
            status = "FRAGMENT"
        rows.append((r.species, r.protein, r.family, r.category, status,
                     len(req), n_ok, round(min(covs.values()) if covs else 0.0, 3),
                     round(sum(covs.values())/len(covs), 3) if covs else 0.0,
                     r.prot_len, round(lr, 2),
                     ";".join(f"{k}:{v:.2f}" for k, v in covs.items())))
    comp = pd.DataFrame(rows, columns=["species","protein","family","category",
        "status","n_required","n_complete","min_domain_cov","mean_domain_cov",
        "prot_len","length_ratio","domain_coverages"])
    if skipped:
        log(f"WARNING: skipped families absent from the map: {sorted(skipped)}")
    comp.to_csv(os.path.join(ddir, "completeness.tsv"), sep="\t", index=False)

    summ = (comp.groupby(["species","family"]).status
              .value_counts().unstack(fill_value=0).reset_index())
    for s in ("COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"):
        if s not in summ:
            summ[s] = 0
    summ["n"] = summ[["COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"]].sum(axis=1)
    summ["pct_complete"] = (summ.COMPLETE / summ.n * 100).round(1)
    summ.to_csv(os.path.join(ddir, "completeness_by_species_family.tsv"),
                sep="\t", index=False)

    byfam = comp.groupby("family").status.value_counts().unstack(fill_value=0)
    byfam["pct_complete"] = (byfam.get("COMPLETE", 0) /
                             byfam.sum(axis=1) * 100).round(1)
    print("\n=== Domain completeness by family ===")
    print(byfam.sort_values("pct_complete").to_string())
    log("domains done")


# =========================================================================== #
CLAN_RE = re.compile(r"clan=(\S+)")


def cmd_cyp(a):
    """Assign CYPs to clans using a labelled reference set.

    Builds one HMM per clan from the reference sequences, then scores every
    CYP against all clan HMMs and takes the best. `margin` is the bitscore gap
    to the runner-up: small margins are assignments you should not trust, and
    the deep clan splits mean a Drosophila-only reference transfers imperfectly
    to Lepidoptera. Check the tree before believing a marginal call."""
    need_pandas()
    for exe in ("mafft", "hmmbuild", "hmmsearch"):
        if not shutil.which(exe):
            sys.exit(f"ERROR: {exe} not found. module load HMMER MAFFT")
    if not os.path.exists(a.refs):
        sys.exit(f"ERROR: reference FASTA not found: {a.refs}")
    cdir = os.path.join(a.out, "cyp"); os.makedirs(os.path.join(cdir, "clan_hmms"), exist_ok=True)

    fam_faa = os.path.join(a.out, "fasta", f"{a.family}.faa")
    if not os.path.exists(fam_faa):
        sys.exit(f"ERROR: {fam_faa} not found. Run `extract` first.")

    by_clan = {}
    for hdr, seq in read_fasta(a.refs):
        m = CLAN_RE.search(hdr)
        if m:
            by_clan.setdefault(m.group(1), []).append((hdr.split()[0], seq))
    if not by_clan:
        sys.exit(f"ERROR: no 'clan=' tags in {a.refs} headers")
    log("reference clans: " + ", ".join(f"{k}({len(v)})" for k, v in sorted(by_clan.items())))

    hmms = []
    for clan, seqs in sorted(by_clan.items()):
        if len(seqs) < 2:
            log(f"skipping clan {clan}: only {len(seqs)} reference sequence"); continue
        base = os.path.join(cdir, "clan_hmms", clan)
        with open(base + ".faa", "w") as fh:
            for n, s in seqs:
                fh.write(f">{n}\n{s}\n")
        with open(base + ".aln", "w") as fh:
            run(["mafft", "--auto", "--quiet", "--thread", a.threads, base + ".faa"], stdout=fh)
        run(["hmmbuild", "--amino", "-n", clan, base + ".hmm", base + ".aln"],
            stdout=subprocess.DEVNULL)
        hmms.append((clan, base + ".hmm"))

    scores = {}
    for clan, hmm in hmms:
        tbl = os.path.join(cdir, f"hits_{clan}.tbl")
        run(["hmmsearch", "--noali", "-E", "1e-5", "--cpu", a.threads,
             "--tblout", tbl, "-o", os.devnull, hmm, fam_faa])
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                scores.setdefault(f[0], {})[clan] = float(f[5])

    rows = []
    for tgt, sc in scores.items():
        ranked = sorted(sc.items(), key=lambda x: -x[1])
        top, tops = ranked[0]
        margin = tops - ranked[1][1] if len(ranked) > 1 else tops
        sp, prot = tgt.split("|", 1) if "|" in tgt else ("NA", tgt)
        rows.append((sp, prot, top, round(tops, 1), round(margin, 1),
                     "LOW" if margin < 20 else "OK"))
    calls = pd.DataFrame(rows, columns=["species","protein","clan","score",
                                        "margin","confidence"])
    unassigned = set()
    for hdr, _ in read_fasta(fam_faa):
        if hdr.split()[0] not in scores:
            unassigned.add(hdr.split()[0])
    calls.to_csv(os.path.join(cdir, f"{a.family}_clan_calls.tsv"), sep="\t", index=False)

    counts = calls.groupby(["species","clan"]).size().unstack(fill_value=0)
    counts["UNASSIGNED"] = pd.Series(
        {s.split("|")[0]: 1 for s in unassigned}).groupby(level=0).sum() \
        if unassigned else 0
    counts = counts.fillna(0).astype(int)
    counts.to_csv(os.path.join(cdir, f"{a.family}_clan_counts.tsv"), sep="\t")
    print(f"\n=== {a.family} clan counts per species ===")
    print(counts.to_string())
    print(f"\nlow-confidence assignments (margin < 20 bits): "
          f"{(calls.confidence == 'LOW').sum()} of {len(calls)}")
    if unassigned:
        log(f"{len(unassigned)} sequences matched no clan HMM above E=1e-5")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fdir = os.path.join(a.out, "figures"); os.makedirs(fdir, exist_ok=True)
        cl = counts.drop(columns=[c for c in ["UNASSIGNED"] if c in counts])
        frac = cl.div(cl.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        frac = frac.loc[frac.sort_values(list(frac.columns)[0]).index]
        fig, ax = plt.subplots(figsize=(9, max(6, .24 * len(frac))))
        left = np.zeros(len(frac))
        for c_ in frac.columns:
            ax.barh(range(len(frac)), frac[c_].values, left=left, label=c_)
            left += frac[c_].values
        ax.set_yticks(range(len(frac))); ax.set_yticklabels(frac.index, fontsize=6)
        ax.set_xlabel(f"fraction of {a.family}s"); ax.legend(fontsize=7, ncol=len(frac.columns))
        ax.set_title(f"{a.family} clan composition")
        fig.tight_layout(); fig.savefig(os.path.join(fdir, f"{a.family}_clans.pdf")); plt.close(fig)
        log(f"wrote {fdir}/{a.family}_clans.pdf")
    except ImportError:
        pass
    log("cyp done")


# =========================================================================== #
def cmd_trees(a):
    """MAFFT alignment and FastTree per family."""
    for exe in ("mafft", "FastTree"):
        if not shutil.which(exe) and not shutil.which(exe.lower()):
            sys.exit(f"ERROR: {exe} not found. module load MAFFT FastTree")
    ft = "FastTree" if shutil.which("FastTree") else "fasttree"
    fdir = os.path.join(a.out, "fasta")
    tdir = os.path.join(a.out, "trees"); os.makedirs(tdir, exist_ok=True)
    fams = a.families.split(",") if a.families else \
        [os.path.basename(f)[:-4] for f in sorted(glob.glob(os.path.join(fdir, "*.faa")))]

    for fam in fams:
        faa = os.path.join(fdir, f"{fam}.faa")
        if not os.path.exists(faa):
            log(f"skipping {fam}: no {faa}"); continue
        n = sum(1 for _ in read_fasta(faa))
        if n < 4:
            log(f"skipping {fam}: only {n} sequences"); continue
        if n > a.max_seqs:
            log(f"skipping {fam}: {n} sequences exceeds --max-seqs {a.max_seqs}"); continue
        aln, tre = os.path.join(tdir, f"{fam}.aln"), os.path.join(tdir, f"{fam}.tre")
        if os.path.exists(tre) and not a.force:
            log(f"{fam}: tree exists, skipping"); continue
        # --retree 1 for large families: --auto picks an O(N^2) method that is
        # impractical past a few thousand sequences.
        mode = ["--auto"] if n <= 1000 else ["--retree", "1", "--maxiterate", "0"]
        with open(aln, "w") as fh:
            run(["mafft", *mode, "--quiet", "--anysymbol",
                 "--thread", a.threads, faa], stdout=fh)
        with open(tre, "w") as fh:
            run([ft, "-lg", "-quiet", aln], stdout=fh)
        log(f"{fam}: {n} sequences -> {tre}")
    log("trees done")


# =========================================================================== #
#  Newick parsing and tree layout. Self-contained: no ete3, no ggtree.
# =========================================================================== #
class Node:
    __slots__ = ("name","length","children","parent","x","y","angle","r","_sz")

    def __init__(self, name="", length=0.0):
        self.name, self.length, self.children, self.parent = name, length, [], None
        self.x = self.y = self.angle = self.r = 0.0
        self._sz = 0

    def is_leaf(self):
        return not self.children


def parse_newick(text):
    text = re.sub(r"\[[^\]]*\]", "", text.strip())     # drop comments
    text = text.rstrip(";").strip()
    pos = [0]

    def parse_node():
        node = Node()
        if text[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(parse_node())
                node.children[-1].parent = node
                if text[pos[0]] == ",":
                    pos[0] += 1
                else:
                    break
            assert text[pos[0]] == ")", f"malformed Newick at {pos[0]}"
            pos[0] += 1
        m = re.match(r"([^,():;]*)", text[pos[0]:])
        label = m.group(1).strip()
        pos[0] += m.end()
        if text[pos[0]:pos[0]+1] == ":":
            pos[0] += 1
            m2 = re.match(r"[-\d.eE+]+", text[pos[0]:])
            node.length = float(m2.group(0)); pos[0] += m2.end()
        # A label on an internal node is support, not a name.
        if node.is_leaf():
            node.name = label.strip("'\"")
        return node
    return parse_node()


def tree_nodes(root):
    stack, out = [root], []
    while stack:
        n = stack.pop(); out.append(n); stack.extend(n.children)
    return out


def tree_leaves(root):
    return [n for n in tree_nodes(root) if n.is_leaf()]


def prune_tree(root, keep):
    """Drop tips not in `keep`, then collapse the resulting single-child nodes."""
    def rec(n):
        if n.is_leaf():
            return n if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        n.children = kids
        for k in kids:
            k.parent = n
        return n
    r = rec(root)
    if r is None:
        sys.exit("ERROR: no tips in the tree match the count matrix")
    r.length = 0.0
    return r


def ladderize(root):
    def rec(n):
        n._sz = 1 if n.is_leaf() else sum(rec(c) for c in n.children)
        n.children.sort(key=lambda c: c._sz)
        return n._sz
    rec(root)
    return root


def set_depths(root, cladogram=True):
    """Radial/horizontal coordinate for every node."""
    if cladogram:
        def depth(n, d=0):
            n.r = d
            for c in n.children:
                depth(c, d + 1)
        depth(root)
        maxd = max(n.r for n in tree_nodes(root))
        # push all tips to the rim so labels line up
        for n in tree_nodes(root):
            n.r = n.r / maxd if maxd else 0.0
        for n in tree_leaves(root):
            n.r = 1.0
    else:
        def dist(n, d=0.0):
            n.r = d
            for c in n.children:
                dist(c, d + max(c.length, 0.0))
        dist(root)
        maxd = max(n.r for n in tree_nodes(root)) or 1.0
        for n in tree_nodes(root):
            n.r /= maxd
    return root


def set_angles(root, span=350.0, start=0.0):
    lv = tree_leaves(root)
    step = span / max(len(lv), 1)
    for i, n in enumerate(lv):
        n.angle = math.radians(start + i * step)
    def rec(n):
        if n.is_leaf():
            return n.angle
        a = [rec(c) for c in n.children]
        n.angle = (min(a) + max(a)) / 2
        return n.angle
    rec(root)
    return root, step


# =========================================================================== #
# Distinct palettes so concentric rings never share a colour. Sharing one made
# CYP clan "CYP3" and diet "Polyphagous" the same orange, which is unreadable.
_PALETTES = [
    ["#1b4965","#2a9d8f","#8ecae6","#023047","#5fa8d3","#0b6e4f","#89c2d9"],   # blues/teals
    ["#bc4749","#f4a261","#e9c46a","#6a040f","#dda15e","#9c6644","#ffb703"],   # warm
    ["#7b2cbf","#c77dff","#5a189a","#e0aaff","#3c096c","#9d4edd","#240046"],   # purples
]

def _palette(levels, ring=0):
    base = _PALETTES[ring % len(_PALETTES)]
    return {lv: base[i % len(base)] for i, lv in enumerate(levels)}


def cmd_tree(a):
    """Radial and rectangular species cladograms with trait and defensome rings.

    Reads any Newick (OrthoFinder's SpeciesTree_rooted.txt works directly),
    prunes it to the species in counts.tsv, then draws:
      - tip labels coloured by the metadata column given to --group-by
      - concentric rings, one per defensome family, z-scored across species
      - an outer bar ring for total CORE defensome
    """
    need_pandas()
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge, Rectangle
    except ImportError:
        sys.exit("ERROR: matplotlib required. pip install --user matplotlib")

    counts = pd.read_csv(os.path.join(a.out, "counts.tsv"), sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    use = counts if a.raw or norm is None else norm
    lab = "copies" if use is counts else f"copies per 10k {size_unit(a.out)}"
    dm = read_map(a.map)

    root = ladderize(prune_tree(parse_newick(open(a.tree).read()), set(counts.index)))
    tips = [n.name for n in tree_leaves(root)]
    log(f"tree: {len(tips)} tips after pruning to counts.tsv")
    missing = sorted(set(counts.index) - set(tips))
    if missing:
        log(f"WARNING: not in the tree: {missing}")

    grp, cmap = None, {}
    if a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        want, key, bestn = set(counts.index), None, 0
        for c_ in md.columns:
            v = md[c_].astype(str).str.strip().str.replace(" ", "_", regex=False)
            if len(want & set(v)) > bestn:
                key, bestn = c_, len(want & set(v))
        if key and a.group_by in md.columns:
            md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
            grp = md.set_index(key)[a.group_by]
            cmap = _palette(sorted(grp.dropna().unique()))
            log(f"trait '{a.group_by}' from column '{key}' ({bestn} matched)")

    fams = (a.families.split(",") if a.families
            else [f for f in dm.loc[dm.tier == "CORE", "family"]
                  if f in use.columns and use[f].sum() > 0])
    fams = [f for f in fams if f in use.columns]
    m = use.loc[tips, fams]
    z = ((m - m.mean()) / m.std().replace(0, np.nan)).fillna(0).clip(-2.5, 2.5)
    total = use.loc[tips, fams].sum(axis=1)

    fig_dir = os.path.join(a.out, "figures"); os.makedirs(fig_dir, exist_ok=True)

    # ------------------------------------------------------------------ radial
    set_depths(root, cladogram=not a.phylogram)
    # Leave a wedge open at the top. Ring labels sit in it: at 90 degrees the
    # radial direction is vertical, so successive rings stack downward-to-upward
    # and horizontal label text never overlaps.
    GAP_DEG = 26.0
    root, step = set_angles(root, span=360.0 - GAP_DEG, start=90.0 + GAP_DEG/2)
    R_TREE, RING_W, GAP = 1.0, 0.085, 0.30
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.set_aspect("equal"); ax.axis("off")

    def pol(r, th):
        return r*math.cos(th), r*math.sin(th)

    for n in tree_nodes(root):                       # branches
        for c in n.children:
            r0 = n.r * R_TREE
            x0, y0 = pol(r0, c.angle); x1, y1 = pol(c.r*R_TREE, c.angle)
            ax.plot([x0, x1], [y0, y1], c="#444", lw=.9, zorder=1)
        if n.children:                                # connector arc
            a0 = min(c.angle for c in n.children); a1 = max(c.angle for c in n.children)
            arc = np.linspace(a0, a1, 40)
            ax.plot(n.r*R_TREE*np.cos(arc), n.r*R_TREE*np.sin(arc),
                    c="#444", lw=.9, zorder=1)

    half = math.radians(step) / 2 * 0.92
    cm_ = plt.cm.RdBu_r
    for k, fam in enumerate(fams):                    # family rings
        r0 = R_TREE + GAP + k*RING_W
        for n in tree_leaves(root):
            v = (z.loc[n.name, fam] + 2.5) / 5.0
            ax.add_patch(Wedge((0, 0), r0+RING_W*.9, math.degrees(n.angle-half),
                               math.degrees(n.angle+half), width=RING_W*.9,
                               facecolor=cm_(v), edgecolor="none", zorder=2))
        ax.text(0.0, r0 + RING_W*.45, fam, fontsize=7, ha="center",
                va="center", color="#333",
                bbox=dict(fc="white", ec="none", pad=.6, alpha=.85), zorder=4)

    r_bar = R_TREE + GAP + len(fams)*RING_W + 0.10    # outer total bar ring
    bmax = total.max() or 1
    bmin = total.min()
    ax.add_patch(Wedge((0, 0), r_bar, 0, 360, width=.004,
                       facecolor="#ccc", edgecolor="none", zorder=1))
    ax.text(0.0, r_bar + 0.16, f"total CORE ({lab})", fontsize=7.5,
            ha="center", va="center", color="#333",
            bbox=dict(fc="white", ec="none", pad=.6, alpha=.85), zorder=4)
    for n in tree_leaves(root):
        h = 0.30 * (total[n.name] - bmin*.85) / (bmax - bmin*.85)
        col = cmap.get(grp.get(n.name), "#888") if grp is not None else "#2a6f97"
        ax.add_patch(Wedge((0, 0), r_bar+h, math.degrees(n.angle-half),
                           math.degrees(n.angle+half), width=h,
                           facecolor=col, edgecolor="none", alpha=.85, zorder=2))
    for n in tree_leaves(root):                       # tip labels
        col = cmap.get(grp.get(n.name), "#222") if grp is not None else "#222"
        deg = math.degrees(n.angle) % 360
        flip = 90 < deg < 270
        ax.text(*pol(r_bar + 0.36, n.angle), n.name.replace("_", " "),
                rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                ha="right" if flip else "left", va="center",
                fontsize=7, color=col, style="italic", zorder=3)

    lim = r_bar + 1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    if grp is not None:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=v, label=str(k_)) for k_, v in cmap.items()],
                  loc="upper left", frameon=False, fontsize=10, title=a.group_by)
    sm = plt.cm.ScalarMappable(cmap=cm_, norm=plt.Normalize(-2.5, 2.5))
    cb = fig.colorbar(sm, ax=ax, shrink=.25, pad=.02, location="right")
    cb.set_label(f"z-score, {lab}", fontsize=9)
    ax.set_title(f"Defensome across the species tree  ({lab})", fontsize=14, pad=18)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "10_tree_radial.pdf"))
    plt.close(fig)

    # ------------------------------------------------------------- rectangular
    lv = tree_leaves(root)
    for i, n in enumerate(lv):
        n.y = i
    def sety(n):
        if n.is_leaf():
            return n.y
        ys = [sety(c) for c in n.children]
        n.y = (min(ys) + max(ys)) / 2
        return n.y
    sety(root)
    fig, (axt, axh) = plt.subplots(1, 2, figsize=(17, max(8, .30*len(lv))),
                                   gridspec_kw={"width_ratios": [1, 2.1]})
    for n in tree_nodes(root):
        for c in n.children:
            axt.plot([n.r, c.r], [c.y, c.y], c="#444", lw=.9)
        if n.children:
            axt.plot([n.r, n.r], [min(c.y for c in n.children),
                                  max(c.y for c in n.children)], c="#444", lw=.9)
    for n in lv:
        col = cmap.get(grp.get(n.name), "#222") if grp is not None else "#222"
        axt.text(n.r + .02, n.y, n.name.replace("_", " "), va="center",
                 fontsize=7, color=col, style="italic")
    axt.set_xlim(-.02, 1.55); axt.set_ylim(-1, len(lv)); axt.axis("off")
    axt.set_title("species tree", fontsize=11)
    im = axh.imshow(z.values, aspect="auto", cmap=cm_, vmin=-2.5, vmax=2.5,
                    extent=[0, len(fams), len(lv)-.5, -.5])
    axh.set_xticks(np.arange(len(fams))+.5)
    axh.set_xticklabels(fams, rotation=90, fontsize=7.5)
    axh.set_yticks([]); axh.set_title(f"defensome, z-scored ({lab})", fontsize=11)
    fig.colorbar(im, ax=axh, shrink=.4, label="z-score")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "11_tree_heatmap.pdf"))
    plt.close(fig)
    log(f"wrote {fig_dir}/10_tree_radial.pdf and 11_tree_heatmap.pdf")


# =========================================================================== #
def cmd_genetree(a):
    """Radial gene tree for one family, annotated with clan and species trait.

    A 4000-tip tree with unlabelled dots tells you nothing. This draws two
    concentric rings (clan, then species trait), labels clan blocks that form
    contiguous arcs, and writes species labels when the tree is small enough
    for them to be legible. It also writes a clan-by-trait contingency table,
    which is usually the number you actually wanted from the picture.
    """
    need_pandas()
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge, Patch
    except ImportError:
        sys.exit("ERROR: matplotlib required")
    tre = a.tree or os.path.join(a.out, "trees", f"{a.family}.tre")
    if not os.path.exists(tre):
        sys.exit(f"ERROR: {tre} not found. Run `trees` first.")
    root = ladderize(parse_newick(open(tre).read()))
    lv = tree_leaves(root)
    n = len(lv)
    log(f"{a.family} gene tree: {n} tips")

    sp_of = {x.name: x.name.split("|")[0] for x in lv}

    # clan, if it has been computed
    clan_of, clan_cmap = {}, {}
    clan_f = os.path.join(a.out, "cyp", f"{a.family}_clan_calls.tsv")
    if os.path.exists(clan_f):
        cl = pd.read_csv(clan_f, sep="\t")
        clan_of = {f"{r.species}|{r.protein}": r.clan for r in cl.itertuples()}
        clan_cmap = _palette(sorted(set(clan_of.values())), ring=0)

    # species trait; with a sample sheet and no other trait, colour by method,
    # which is the variable a method-comparison run is actually about
    trait_of, trait_cmap, trait_name = {}, {}, None
    if getattr(a, "samplesheet", None) and not (a.metadata and a.group_by):
        sh_ = pd.read_csv(a.samplesheet, sep="\t")
        mth = dict(zip(sh_.sample_id.astype(str), sh_.method.astype(str)))
        trait_of = {t: mth.get(sp_of[t]) for t in sp_of}
        if any(trait_of.values()):
            trait_cmap = _palette(sorted({v for v in trait_of.values() if v}), ring=1)
            trait_name = "method"
            log(f"colouring by method from {a.samplesheet}")
    elif a.metadata and a.group_by:
        md = pd.read_csv(a.metadata, sep="\t", comment="#")
        sp_all = set(sp_of.values())
        key, bestn = None, 0
        for c_ in md.columns:
            v = md[c_].astype(str).str.strip().str.replace(" ", "_", regex=False)
            if len(sp_all & set(v)) > bestn:
                key, bestn = c_, len(sp_all & set(v))
        if key and a.group_by in md.columns:
            md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
            g = md.set_index(key)[a.group_by]
            trait_of = {t: g.get(sp_of[t]) for t in sp_of}
            trait_cmap = _palette(sorted({v for v in trait_of.values() if pd.notna(v)}), ring=1)
            trait_name = a.group_by
            log(f"trait '{a.group_by}' from column '{key}' ({bestn} species matched)")

    rings = []
    if clan_of:
        rings.append(("clan", clan_of, clan_cmap))
    if trait_of:
        rings.append((trait_name, trait_of, trait_cmap))
    if not rings:                      # fall back to colouring by species
        cm_ = _palette(sorted(set(sp_of.values())), ring=2)
        rings.append(("species", sp_of, cm_))

    set_depths(root, cladogram=False)
    GAP_DEG = 16.0
    root, step = set_angles(root, span=360.0 - GAP_DEG, start=90.0 + GAP_DEG/2)

    R, RW = 1.0, 0.055
    fig, ax = plt.subplots(figsize=(16, 16)); ax.set_aspect("equal"); ax.axis("off")
    for nd in tree_nodes(root):
        for c in nd.children:
            ax.plot([nd.r*R*math.cos(c.angle), c.r*R*math.cos(c.angle)],
                    [nd.r*R*math.sin(c.angle), c.r*R*math.sin(c.angle)],
                    c="#666", lw=.3, zorder=1)
        if nd.children:
            arc = np.linspace(min(c.angle for c in nd.children),
                              max(c.angle for c in nd.children), 24)
            ax.plot(nd.r*R*np.cos(arc), nd.r*R*np.sin(arc), c="#666", lw=.3, zorder=1)

    half = math.radians(step)/2 * 1.0
    for k, (name, mapping, cm_) in enumerate(rings):
        r0 = R + 0.06 + k*RW
        for nd in lv:
            col = cm_.get(mapping.get(nd.name), "#e9ecef")
            ax.add_patch(Wedge((0, 0), r0+RW*.92, math.degrees(nd.angle-half),
                               math.degrees(nd.angle+half), width=RW*.92,
                               facecolor=col, edgecolor="none", zorder=2))
        ax.text(0.0, r0+RW*.45, name, fontsize=9, ha="center", va="center",
                bbox=dict(fc="white", ec="#ddd", pad=.9, alpha=.95,
                          boxstyle="round,pad=0.3"), zorder=5)

    r_out = R + 0.06 + len(rings)*RW

    # label clan blocks that form contiguous arcs
    if clan_of:
        seq = [clan_of.get(nd.name) for nd in lv]
        i = 0
        min_run = max(12, n // 60)
        while i < len(seq):
            j = i
            while j+1 < len(seq) and seq[j+1] == seq[i]:
                j += 1
            if seq[i] is not None and (j - i + 1) >= min_run:
                a0, a1 = lv[i].angle, lv[j].angle
                mid = (a0 + a1) / 2
                arc = np.linspace(a0-half, a1+half, 40)
                rr = r_out + 0.05
                ax.plot(rr*np.cos(arc), rr*np.sin(arc),
                        c=clan_cmap.get(seq[i], "#333"), lw=3.0, solid_capstyle="butt",
                        zorder=3)
                deg = math.degrees(mid) % 360
                flip = 90 < deg < 270
                ax.text((rr+0.05)*math.cos(mid), (rr+0.05)*math.sin(mid),
                        f"{seq[i]} ({j-i+1})",
                        rotation=deg+180 if flip else deg, rotation_mode="anchor",
                        ha="right" if flip else "left", va="center", fontsize=8,
                        color=clan_cmap.get(seq[i], "#333"), zorder=4)
            i = j + 1
        r_out += 0.30

    # species labels only when they can actually be read
    label_limit = a.max_labels
    if n <= label_limit:
        for nd in lv:
            deg = math.degrees(nd.angle) % 360
            flip = 90 < deg < 270
            col = trait_cmap.get(trait_of.get(nd.name), "#333") if trait_of else "#333"
            ax.text((r_out+0.02)*math.cos(nd.angle), (r_out+0.02)*math.sin(nd.angle),
                    sp_of[nd.name].replace("_", " "),
                    rotation=deg+180 if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center",
                    fontsize=max(2.0, min(6.0, 260.0/n)), color=col, style="italic",
                    zorder=4)
        r_out += 0.42
    else:
        ax.text(0, -r_out-0.16,
                f"{n} tips: species labels suppressed above --max-labels {label_limit}. "
                f"Ring colours and {a.family}_tip_table.tsv carry the same information.",
                ha="center", fontsize=9, color="#666")
        r_out += 0.22

    handles = []
    for name, mapping, cm_ in rings:
        handles.append(Patch(color="white", label=f"— {name} —"))
        handles += [Patch(color=v, label=str(k_)) for k_, v in cm_.items()][:12]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5,
              bbox_to_anchor=(-0.02, 1.0))
    ax.set_xlim(-r_out-.12, r_out+.12); ax.set_ylim(-r_out-.12, r_out+.12)
    ttl = f"{a.family} gene tree: {n} sequences from {len(set(sp_of.values()))} species"
    ax.set_title(ttl, fontsize=15, pad=16)
    fig_dir = os.path.join(a.out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    tag = "_".join(r[0] for r in rings).replace(" ", "")
    out = os.path.join(fig_dir, f"12_{a.family}_genetree_{tag}.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    log(f"wrote {out}")

    # the numbers behind the picture
    tips = pd.DataFrame({
        "tip": [nd.name for nd in lv],
        "species": [sp_of[nd.name] for nd in lv],
        "clan": [clan_of.get(nd.name) for nd in lv],
        "trait": [trait_of.get(nd.name) for nd in lv],
        "tip_order": range(len(lv)),
    })
    tips.to_csv(os.path.join(a.out, "trees", f"{a.family}_tip_table.tsv"),
                sep="\t", index=False)
    if clan_of:
        ct = pd.crosstab(tips.species, tips.clan)
        ct["total"] = ct.sum(axis=1)
        ct.to_csv(os.path.join(a.out, "trees", f"{a.family}_clan_by_species.tsv"), sep="\t")
        print(f"\n=== {a.family} clan counts per species ===")
        print(ct.to_string())
        if trait_of:
            ct2 = pd.crosstab(tips.trait, tips.clan)
            print(f"\n=== {a.family} clan by {trait_name} (totals across species) ===")
            print(ct2.to_string())
            frac = ct2.div(ct2.sum(axis=1), axis=0).round(3)
            print("\nas fractions within each group:"); print(frac.to_string())
            ct2.to_csv(os.path.join(a.out, "trees", f"{a.family}_clan_by_{trait_name}.tsv"),
                       sep="\t")

            # per-species clan composition, which is the comparable version
            per = pd.crosstab(tips.species, tips.clan)
            sz = os.path.join(a.out, "proteome_sizes.tsv")
            if os.path.exists(sz):
                nprot = pd.read_csv(sz, sep="\t").set_index("species").n_proteins
                per10k = per.div(nprot.reindex(per.index), axis=0) * 10000
                per10k.round(2).to_csv(
                    os.path.join(a.out, "trees", f"{a.family}_clan_per10k.tsv"), sep="\t")
                tmap = tips.drop_duplicates("species").set_index("species").trait
                fig, ax = plt.subplots(figsize=(10, max(6, .26*len(per10k))))
                order = per10k.sum(axis=1).sort_values().index
                left = np.zeros(len(order))
                for c_ in per10k.columns:
                    ax.barh(range(len(order)), per10k.loc[order, c_].values, left=left,
                            label=c_, color=clan_cmap.get(c_, "#999"))
                    left += per10k.loc[order, c_].values
                ax.set_yticks(range(len(order)))
                ax.set_yticklabels([s.replace("_", " ") for s in order], fontsize=7.5,
                                   style="italic")
                if trait_cmap:
                    for t_, s_ in zip(ax.get_yticklabels(), order):
                        t_.set_color(trait_cmap.get(tmap.get(s_), "#222"))
                ax.set_xlabel(f"{a.family} copies per 10k {size_unit(a.out)}")
                ax.legend(fontsize=8, ncol=len(per10k.columns), loc="lower right")
                ax.set_title(f"{a.family} clan composition per species "
                             f"(tip labels coloured by {trait_name})")
                fig.tight_layout()
                fig.savefig(os.path.join(fig_dir, f"13_{a.family}_clans_per_species.pdf"))
                plt.close(fig)
                log(f"wrote {fig_dir}/13_{a.family}_clans_per_species.pdf")


# =========================================================================== #
def size_unit(out):
    """What the per-10k denominator actually counts.

    After `collapse` the proteomes hold one sequence per gene, so the
    denominator is genes. Without that step it is proteins, and with several
    isoforms per gene those are not the same number. Labelling both "proteins"
    overstated what the ratio means."""
    return "genes" if os.path.exists(os.path.join(out, "collapse_stats.tsv")) else "proteins"


def _read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw) if os.path.exists(path) else None


def build_payload(out, mapfile, metadata=None, group_by=None, full=True,
                  samplesheet=None):
    """Collect every output file that exists into one JSON-serialisable dict."""
    dm = read_map(mapfile)
    counts = _read_tsv(os.path.join(out, "counts.tsv"), index_col=0)
    if counts is None:
        sys.exit(f"ERROR: {out}/counts.tsv not found. Run `annotate` first.")
    per10k = _read_tsv(os.path.join(out, "counts_per10k.tsv"), index_col=0)
    sizes = _read_tsv(os.path.join(out, "proteome_sizes.tsv"))
    qc = _read_tsv(os.path.join(out, "qc_annotation.tsv"), index_col=0)
    qcf = _read_tsv(os.path.join(out, "qc_families.tsv"), index_col=0)
    qcz = _read_tsv(os.path.join(out, "qc_zero_families.tsv"))
    fsum = _read_tsv(os.path.join(out, "family_summary.tsv"))
    comp = _read_tsv(os.path.join(out, "domains", "completeness.tsv"))
    arch = _read_tsv(os.path.join(out, "domains", "architectures.tsv"))
    calls = _read_tsv(os.path.join(out, "gene_calls_all.tsv"))

    species = list(counts.index)
    fams = [f for f in dm.family if f in counts.columns]

    trait, trait_name = {}, None
    if metadata and group_by and os.path.exists(metadata):
        md = pd.read_csv(metadata, sep="\t", comment="#")
        if group_by in md.columns:
            key, best = None, 0
            want = set(species)
            for c in md.columns:
                v = md[c].astype(str).str.strip().str.replace(" ", "_", regex=False)
                if len(want & set(v)) > best:
                    key, best = c, len(want & set(v))
            if key:
                md[key] = md[key].astype(str).str.strip().str.replace(" ", "_", regex=False)
                g = md.set_index(key)[group_by]
                trait = {s: (None if pd.isna(g.get(s)) else str(g.get(s))) for s in species}
                trait_name = group_by
                log(f"dashboard trait '{group_by}' from column '{key}' ({best} matched)")

    samples = {}
    if samplesheet and os.path.exists(samplesheet):
        sh_ = pd.read_csv(samplesheet, sep="\t")
        for r in sh_.itertuples(index=False):
            d = r._asdict()
            samples[str(d["sample_id"])] = {
                "species": str(d.get("species", "") or ""),
                "method": str(d.get("method", "") or ""),
                "pipeline": str(d.get("pipeline", "") or "")}
        if not trait and any(v["method"] for v in samples.values()):
            trait = {s_: samples.get(s_, {}).get("method") or None for s_ in species}
            trait_name = "method"
            log("dashboard trait defaults to 'method' from the sample sheet")

    nprot = {}
    if sizes is not None:
        nprot = dict(zip(sizes.species, sizes.n_proteins.astype(int)))

    sp_rows = []
    for s in species:
        row = {"name": s, "n_proteins": int(nprot.get(s, 0)),
               "trait": trait.get(s),
               "species_name": samples.get(s, {}).get("species", ""),
               "method": samples.get(s, {}).get("method", ""),
               "pipeline": samples.get(s, {}).get("pipeline", "")}
        if qc is not None and s in qc.index:
            row["max_ratio"] = float(qc.loc[s].get("max_ratio", float("nan")))
            row["flag"] = str(qc.loc[s].get("flag", ""))
            row["worst_family"] = str(qc.loc[s].get("worst_family", ""))
        sp_rows.append(row)

    P = {
        "size_unit": size_unit(out),
        "generated": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "out_dir": os.path.abspath(out),
        "trait_name": trait_name,
        "species": sp_rows,
        "families": [{"family": r.family, "category": r.category, "tier": r.tier,
                      "pfam": ",".join(r.pfam_ids), "rule": r.rule,
                      "min_cov": float(r.min_cov), "min_len": int(r.min_len),
                      "notes": ("" if pd.isna(r.notes) else str(r.notes))}
                     for r in dm.itertuples() if r.family in counts.columns],
        "counts": {f: [int(v) for v in counts[f]] for f in fams},
        "per10k": ({f: [round(float(v), 4) for v in per10k[f]] for f in fams}
                   if per10k is not None else None),
    }

    if fsum is not None:
        P["family_summary"] = fsum.where(pd.notna(fsum), None).to_dict("records")
    if qcf is not None:
        P["qc_families"] = [dict(family=i, **{k: (None if pd.isna(v) else v)
                                              for k, v in r.items()})
                            for i, r in qcf.iterrows()]
    if qcz is not None:
        P["qc_zero"] = qcz.where(pd.notna(qcz), None).to_dict("records")
    cr = _read_tsv(os.path.join(out, "counts_rescued.tsv"), index_col=0)
    rcalls = _read_tsv(os.path.join(out, "rescue_calls.tsv"))
    P["rescued"] = {}
    if cr is not None and rcalls is not None and len(rcalls):
        for f in sorted(rcalls[rcalls.validation == "PASS"].family.unique()):
            if f in cr.columns and f in counts.columns:
                extra = (cr[f] - counts[f]).reindex(counts.index).fillna(0)
                P["rescued"][f] = [int(v) for v in extra]

    if comp is not None:
        t = comp.groupby(["family", "status"]).size().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in t:
                t[s_] = 0
        P["completeness"] = {f: {k: int(t.loc[f, k]) for k in
                                 ["COMPLETE","PARTIAL","FRAGMENT","MISSING_DOMAIN"]}
                             for f in t.index}
        sf = comp.groupby(["species", "family"]).status.value_counts().unstack(fill_value=0)
        for s_ in ("COMPLETE", "PARTIAL", "FRAGMENT", "MISSING_DOMAIN"):
            if s_ not in sf:
                sf[s_] = 0
        P["comp_sf"] = [[i[0], i[1], int(r.COMPLETE), int(r.PARTIAL),
                         int(r.FRAGMENT), int(r.MISSING_DOMAIN)]
                        for i, r in sf.iterrows()]
        # the non-complete proteins are the ones worth eyeballing
        bad = comp[comp.status != "COMPLETE"]
        P["flagged_proteins"] = bad[["species","protein","family","status","prot_len",
                                     "length_ratio","min_domain_cov","domain_coverages"]] \
            .head(4000).values.tolist()

    if arch is not None and calls is not None:
        a2 = arch.merge(calls[["species","protein","family"]].drop_duplicates(),
                        on=["species","protein"], how="left")
        ag = (a2.groupby(["family","architecture"]).size()
                .reset_index(name="n").sort_values(["family","n"], ascending=[True, False]))
        P["architectures"] = ag.values.tolist()

    P["stats"] = {}
    for f in glob.glob(os.path.join(out, "kruskal_*.tsv")):
        t = pd.read_csv(f, sep="\t")
        P["stats"][os.path.basename(f)[8:-4]] = t.where(pd.notna(t), None).to_dict("records")
    P["group_means"] = {}
    for f in glob.glob(os.path.join(out, "group_means_*.tsv")):
        t = pd.read_csv(f, sep="\t", index_col=0)
        P["group_means"][os.path.basename(f)[12:-4]] = {
            "groups": [c for c in t.columns if c != "n_groups"],
            "data": {i: [None if pd.isna(v) else float(v)
                         for k, v in r.items() if k != "n_groups"]
                     for i, r in t.iterrows()}}

    # provenance: what was found, what was not, so a reader can see the gaps
    prov = []
    for label, rel in [("hmmsearch domtblout", "hmmsearch"),
                       ("gene calls", "gene_calls_all.tsv"),
                       ("count matrix", "counts.tsv"),
                       ("normalised counts", "counts_per10k.tsv"),
                       ("proteome sizes", "proteome_sizes.tsv"),
                       ("annotation QC", "qc_annotation.tsv"),
                       ("domain table", "domains/domain_table.tsv"),
                       ("completeness", "domains/completeness.tsv"),
                       ("architectures", "domains/architectures.tsv"),
                       ("per-family FASTA", "fasta"),
                       ("gene trees", "trees"),
                       ("CYP clan calls", "cyp"),
                       ("figures", "figures")]:
        p_ = os.path.join(out, rel)
        found = os.path.exists(p_)
        n = len(glob.glob(os.path.join(p_, "*"))) if found and os.path.isdir(p_) else None
        prov.append({"item": label, "path": rel, "found": bool(found),
                     "n_files": n,
                     "modified": (__import__("datetime").datetime.fromtimestamp(
                         os.path.getmtime(p_)).isoformat(timespec="seconds")
                         if found else None)})
    P["provenance"] = prov

    clan_dir = os.path.join(out, "trees")
    P["clans"] = {}
    for f in glob.glob(os.path.join(clan_dir, "*_clan_by_species.tsv")):
        fam = os.path.basename(f).replace("_clan_by_species.tsv", "")
        t = pd.read_csv(f, sep="\t", index_col=0)
        t = t.drop(columns=[c for c in ["total"] if c in t.columns])
        P["clans"][fam] = {"clans": list(t.columns),
                           "data": {i: [int(v) for v in r] for i, r in t.iterrows()}}

    P["samples"] = samples
    cmpd = os.path.join(out, "comparison")
    def _recs(p, **kw):
        t = _read_tsv(p, **kw)
        return None if t is None else t.where(pd.notna(t), None).to_dict("records")
    P["cmp"] = {
        "long": _recs(os.path.join(cmpd, "long_counts.tsv")),
        "recovery": _recs(os.path.join(cmpd, "recovery.tsv")),
        "tests": _recs(os.path.join(cmpd, "paired_tests.tsv")),
        "collapse": _recs(os.path.join(out, "collapse_stats.tsv")),
    }
    so = _read_tsv(os.path.join(out, "short_orf_summary.tsv"))
    P["short_orfs"] = None if so is None else so.to_dict("records")
    cc = _read_tsv(os.path.join(out, "composition_candidates.tsv"))
    P["mt_composition"] = ({} if cc is None or not len(cc)
                           else cc.groupby("species").size().astype(int).to_dict())
    P["tip_clan"] = {}
    for f in glob.glob(os.path.join(out, "cyp", "*_clan_calls.tsv")):
        t = pd.read_csv(f, sep="\t")
        for r in t.itertuples(index=False):
            P["tip_clan"][f"{r.species}|{r.protein}"] = str(r.clan)
    P["tip_status"] = {}
    if comp is not None:
        for r in comp[comp.status != "COMPLETE"].itertuples(index=False):
            P["tip_status"][f"{r.species}|{r.protein}"] = str(r.status)

    P["trees"] = {}
    for label, path in [("species", os.path.join(out, "species_tree.nwk"))]:
        if os.path.exists(path):
            P["trees"][label] = open(path).read().strip()
    for f in sorted(glob.glob(os.path.join(out, "trees", "*.tre"))):
        nwk = open(f).read().strip()
        if len(nwk) < 4_000_000:
            P["trees"][os.path.basename(f)[:-4]] = nwk

    if full and calls is not None:
        fidx = {f: i for i, f in enumerate(fams)}
        sidx = {s: i for i, s in enumerate(species)}
        P["gene_calls"] = [[sidx[r.species], fidx.get(r.family, -1), r.protein,
                            int(r.prot_len), round(float(r.hmm_cov), 3)]
                           for r in calls.itertuples()
                           if r.species in sidx and r.family in fidx]
    return P


def cmd_dashboard(a):
    """Single self-contained HTML dashboard. No server, no CDN, no internet."""
    need_pandas()
    import json
    P = build_payload(a.out, a.map, a.metadata, a.group_by, full=not a.light,
                      samplesheet=a.samplesheet)
    if not P["tip_clan"] and os.path.isdir(os.path.join(a.out, "trees")):
        log("NOTE: no CYP clan calls found, so the tree clan ring will be empty. "
            "Run `cyp --refs <clan-labelled FASTA>` first.")
    if a.species_tree and os.path.exists(a.species_tree):
        P["trees"]["species"] = open(a.species_tree).read().strip()
    P["tool_version"] = __version__
    html = (asset("dashboard.html")
            .replace("/*__JS__*/", asset("dashboard.js"))
            .replace("/*__PAYLOAD__*/null",
                     json.dumps(P, separators=(",", ":"), default=str)))
    dest = a.dashboard_out or os.path.join(a.out, "dashboard.html")
    with open(dest, "w") as fh:
        fh.write(html)
    mb = os.path.getsize(dest) / 1e6
    log(f"wrote {dest} ({mb:.1f} MB)")
    print(f"\nOpen it directly in a browser. Nothing else is needed:\n  {os.path.abspath(dest)}")
    if mb > 40:
        print("\nLarge file. Rerun with --light to drop the per-protein tables.")


# =========================================================================== #
#  CYP clan reference building and external-dataset triage
# =========================================================================== #
def norm_clan(v):
    if v is None:
        return None
    s = str(v).strip().upper().replace("CLAN", "").replace("_", "").replace("-", "").strip()
    if s.startswith("MITO"):
        return "MITO"
    m = re.search(r"([234])", s)
    return {"2": "CYP2", "3": "CYP3", "4": "CYP4"}.get(m.group(1)) if m else None


def clan_from_name(name):
    """Fallback: infer the clade from the CYP family number using Feyereisen's
    assignment of insect families to the four clades."""
    m = re.search(r"CYP(\d+)", str(name).upper())
    if not m:
        return None
    f = int(m.group(1))
    if f in (12, 301, 302, 306, 314, 315):
        return "MITO"
    if f in (4, 311, 312, 313, 316, 318, 325, 340, 341, 349, 351, 353):
        return "CYP4"
    if f in (6, 9, 28, 308, 309, 317, 321, 324, 327, 332, 337, 345, 346, 347):
        return "CYP3"
    if f in (15, 18, 303, 304, 305, 307) or 300 <= f <= 307:
        return "CYP2"
    return None


def build_icpd_refs(fasta, table, orders, species, min_len, max_per_clan, out,
                    id_col=None, include_fragments=False):
    """Clan-labelled reference FASTA from the ICPD protein set.

    Uses a clan column from the info table when one exists, and otherwise
    infers the clade from the CYP family number. Caps sequences per clan and
    samples evenly across species, because ten near-identical CYP6 paralogues
    from one over-studied pest would otherwise pull every query toward it."""
    meta = {}
    fasta_ids = {h.split()[0] for h, _ in read_fasta(fasta)} if fasta else set()
    if table:
        need_pandas()
        t = (pd.read_excel(table) if str(table).endswith((".xlsx", ".xls"))
             else pd.read_csv(table, sep=None, engine="python"))
        print(f"info table: {t.shape[0]} rows; columns: {list(t.columns)}")
        # Pick the ID column by testing which one actually matches the FASTA
        # headers. ICPD's FASTA is keyed on "Database ID" (Abter001) while the
        # first ID-looking column is "Protein ID"
        # (Abscondita_terminalisg020.t1); guessing from the column name picks
        # the wrong one and every join silently fails.
        idc, best = None, 0
        if fasta_ids:
            for c in t.columns:
                n = len(fasta_ids & set(t[c].astype(str).str.strip()))
                if n > best:
                    idc, best = c, n
            if idc:
                print(f"id column '{idc}' matches {best}/{len(fasta_ids)} FASTA headers")
        if not idc:
            idc = next((c for c in t.columns if re.search(
                r"^(id|gene.?id|protein.?id|accession|seq.?id)", str(c), re.I)), t.columns[0])
            if fasta_ids:
                print(f"WARNING: no column matched the FASTA headers; falling back to '{idc}'.\n"
                      "         Check --id-col, or drop --fasta and build from the table's "
                      "protein sequence column.")
        spc = next((c for c in t.columns if re.search(r"species|organism", str(c), re.I)), None)
        ordc = next((c for c in t.columns if re.search(r"order", str(c), re.I)), None)
        clc = next((c for c in t.columns if re.search(r"clan|clade", str(c), re.I)), None)
        nmc = next((c for c in t.columns if re.search(r"genesymbol|description", str(c), re.I)),
                   None) or next((c for c in t.columns if re.search(r"symbol|cyp", str(c), re.I)), None)
        cmp_ = next((c for c in t.columns if re.search(r"complete", str(c), re.I)), None)
        seqc = next((c for c in t.columns if re.search(r"protein.*seq", str(c), re.I)), None)
        if id_col:
            idc = id_col
        print(f"using id='{idc}' species='{spc}' order='{ordc}' clan='{clc}' "
              f"name='{nmc}' completeness='{cmp_}' sequence='{seqc}'")
        for r in t.itertuples(index=False):
            d = dict(zip(t.columns, r))
            meta[str(d[idc]).strip()] = {
                "species": str(d.get(spc, "")).strip() if spc else "",
                "order": str(d.get(ordc, "")).strip() if ordc else "",
                "clan": norm_clan(d.get(clc)) if clc else None,
                "name": str(d.get(nmc, "")).strip() if nmc else "",
                "complete": str(d.get(cmp_, "")).strip().lower() if cmp_ else "",
                "seq": str(d.get(seqc, "")).strip() if seqc else ""}
        if cmp_:
            vc = t[cmp_].astype(str).str.strip().str.lower().value_counts()
            print(f"\ncompleteness in the table: {vc.to_dict()}")

    keep_ord = {o.strip().lower() for o in (orders or "").split(",") if o.strip()}
    keep_sp = {sp.strip().lower() for sp in (species or "").split(",") if sp.strip()}
    by_clan, no_clan, filtered, frag, bad_seq = {}, 0, 0, 0, 0
    if fasta:
        src = ((h.split()[0], s_) for h, s_ in read_fasta(fasta))
    else:
        src = ((k, v.get("seq", "")) for k, v in meta.items())
        print("no --fasta given: building from the table's protein sequence column")
    for sid, seq in src:
        m = meta.get(sid) or {"species": sid, "order": sid, "clan": None, "name": sid}
        # ICPD marks a large fraction of its predictions "fragmented". A
        # truncated P450 in a clan reference degrades the profile it goes into,
        # and these are exactly the models that lack the diagnostic motifs.
        # ICPD marks entries "full-length", "fragmented" or "artificially
        # fused". Fragments lack the diagnostic motifs; fused models are two
        # adjacent genes merged into one chimera, which is worse. Neither
        # belongs in a profile HMM.
        cmpl = m.get("complete", "")
        if not include_fragments and (cmpl.startswith("frag") or "fused" in cmpl):
            frag += 1; continue
        # ---- BUG 3: some ICPD "full-length" sequences carry internal stop
        # codons. cd-hit silently discards them; hmmbuild would keep them and
        # corrupt the alignment column statistics.
        core = seq.rstrip("*")
        if "*" in core or "X" * 10 in core.upper():
            bad_seq += 1; continue
        if keep_ord and not any(o in str(m.get("order", "")).lower()
                                or o in str(m.get("species", "")).lower() for o in keep_ord):
            filtered += 1; continue
        if keep_sp and not any(x in str(m.get("species", "")).lower() for x in keep_sp):
            filtered += 1; continue
        if len(seq) < min_len:
            filtered += 1; continue
        hm = CLAN_RE.search(sid)
        clan = (m.get("clan") or (norm_clan(hm.group(0)) if hm else None)
                or clan_from_name(m.get("name") or m.get("description") or sid))
        if not clan:
            no_clan += 1; continue
        by_clan.setdefault(clan, []).append((sid, m.get("species", ""), m.get("name", ""), seq))

    if not by_clan:
        sys.exit("ERROR: nothing could be assigned to a clan. Check the column "
                 "detection printed above, or supply a table with a clan column.")
    d = os.path.dirname(os.path.abspath(out))
    if d:
        os.makedirs(d, exist_ok=True)
    n = 0
    with open(out, "w") as fh:
        for clan, items in sorted(by_clan.items()):
            items.sort(key=lambda x: (x[1], x[0]))
            step = max(1, len(items) // max(max_per_clan, 1))
            for sid, sp_, nm, seq in items[::step][:max_per_clan]:
                fh.write(f">{sid} clan={clan} species={sp_.replace(' ', '_')} name={nm}\n"
                         f"{seq.rstrip('*')}\n")
                n += 1
    print("\n=== reference composition ===")
    for clan, items in sorted(by_clan.items()):
        print(f"  {clan:6s} available {len(items):6d}  written {min(len(items), max_per_clan):5d}")
    print(f"\nskipped: {filtered} filtered out, {frag} fragmented or fused, "
          f"{bad_seq} with internal stops or long X runs, "
          f"{no_clan} with no assignable clan")
    log(f"wrote {n} sequences to {out}")
    print("\nNext, reduce redundancy so one over-sampled pest cannot dominate:")
    print(f"  cd-hit -i {out} -o {out}.c80.faa -c 0.8 -n 5")
    print(f"  python3 defensome.py cyp --out <results> --refs {out}.c80.faa")


def cmd_cyp_refs(a):
    """Build a clan-labelled CYP reference FASTA for `cyp --refs`.

    --source flybase  uses the Drosophila set you already have. Fast, but
                      ~350 Myr from Lepidoptera, so expect LOW-confidence calls.
    --source icpd     uses the Insect Cytochrome P450 Database (66,477 P450s
                      from 682 insect species). Restrict with --orders
                      Lepidoptera and the clan calls get much better.
    """
    if a.source == "flybase":
        if not a.fasta or not os.path.exists(a.fasta):
            sys.exit("ERROR: --fasta is required for --source flybase\n"
                     "  e.g. db/flybase/flybase_dmel_clan_ref.canonical.faa")
        n = 0
        d = os.path.dirname(os.path.abspath(a.out))
        if d:
            os.makedirs(d, exist_ok=True)
        seen = {}
        with open(a.out, "w") as fh:
            for hdr, seq in read_fasta(a.fasta):
                m = re.search(r"clan=(\S+)", hdr)
                clan = norm_clan(m.group(1)) if m else clan_from_name(hdr)
                if not clan or len(seq) < a.min_len:
                    continue
                seen[clan] = seen.get(clan, 0) + 1
                fh.write(f">{hdr.split()[0]} clan={clan} species=Drosophila_melanogaster\n{seq}\n")
                n += 1
        if not n:
            sys.exit(f"ERROR: no clan= tags or recognisable CYP names in {a.fasta}")
        print("=== reference composition ===")
        for k, v in sorted(seen.items()):
            print(f"  {k:6s} {v:5d}")
        log(f"wrote {n} sequences to {a.out}")
        print("\nThis is a Drosophila-only reference. Marginal clan calls will be "
              "flagged LOW.\nFor Lepidoptera, download ICPD and rerun with "
              "--source icpd (see `dataset-help`).")
        return

    if a.fasta and not os.path.exists(a.fasta):
        sys.exit(f"ERROR: --fasta not found: {a.fasta}\n"
                 "  python3 defensome.py dataset-help")
    if not a.fasta and not a.table:
        sys.exit("ERROR: give --fasta, --table, or both.\n"
                 "  python3 defensome.py dataset-help")
    build_icpd_refs(a.fasta, a.table, a.orders, a.species, a.min_len,
                    a.max_per_clan, a.out, a.id_col, a.include_fragments)


INSECT_HINTS = ["insect", "arthropod", "drosophila", "bombyx", "helicoverpa",
                "spodoptera", "plutella", "anopheles", "aedes", "culex",
                "tribolium", "apis", "myzus", "nilaparvata", "papilio",
                "musca", "leptinotarsa", "bemisia", "tetranychus", "acyrthosiphon",
                "animal", "lepidoptera", "diptera", "coleoptera", "hymenoptera"]


def triage_p450rdb(indir, out):
    os.makedirs(out, exist_ok=True)
    def find(pat):
        for f in sorted(os.listdir(indir)):
            if re.search(pat, f, re.I):
                return os.path.join(indir, f)
        return None
    fr = find(r"Reactions")
    if not fr:
        sys.exit(f"ERROR: no Reactions CSV in {indir}")
    rx = pd.read_csv(fr, low_memory=False)
    print(f"Reactions: {len(rx):,} rows x {rx.shape[1]} columns")
    if "Species" in rx.columns:
        k = rx.Species.fillna("unknown").value_counts()
        print("\n=== Rows by kingdom ==="); print(k.to_string())
        k.to_csv(os.path.join(out, "rows_by_kingdom.tsv"), sep="\t")
    sp = rx["Species name"].fillna("unknown").value_counts()
    print(f"\n=== Distinct source organisms: {len(sp)} ==="); print(sp.head(15).to_string())
    sp.to_csv(os.path.join(out, "rows_by_species.tsv"), sep="\t")
    blob = (rx["Species name"].fillna("") + " " +
            rx.get("Species", pd.Series([""] * len(rx))).fillna("")).str.lower()
    ins = rx[blob.apply(lambda s: any(h in s for h in INSECT_HINTS))]
    print(f"\n=== Rows matching an insect/arthropod hint: {len(ins):,} of {len(rx):,} "
          f"({len(ins)/max(len(rx),1)*100:.1f}%) ===")
    if len(ins):
        print(ins["Species name"].value_counts().head(12).to_string())
        ins.to_csv(os.path.join(out, "insect_rows.tsv"), sep="\t", index=False)
    else:
        print("NONE. P450RDB cannot supply insect enzyme-substrate pairs.\n"
              "Use it for host-plant chemistry; get insect P450 data from ICPD.")
    m = melt_reactions(rx)
    m = m[m.sub_SMILES.notna() & ~m.sub_SMILES.astype(str).str.strip().isin(["", "/"])]
    if "sequence" in m.columns:
        m = m[m.sequence.notna() & (m.sequence.astype(str).str.len() > 100)]
    m.to_csv(os.path.join(out, "enzyme_substrate_pairs.tsv"), sep="\t", index=False)
    print(f"\n=== Pairs with both a sequence and a SMILES: {len(m):,} ===")
    print(f"    distinct enzymes: {m['Uniprot ID'].nunique() if 'Uniprot ID' in m else '?'}"
          f"  distinct substrates: {m.substrate.nunique()}")
    print("    POSITIVES ONLY. There are no curated non-substrates here, and")
    print("    absence of a record is not evidence of a non-substrate.")
    if "Transformations" in rx.columns:
        t = rx.Transformations.fillna("unknown").value_counts()
        print("\n=== Transformation types ==="); print(t.head(12).to_string())
        t.to_csv(os.path.join(out, "transformation_types.tsv"), sep="\t")
    if "Species" in rx.columns:
        pl = melt_reactions(rx[rx.Species.fillna("").str.lower().eq("plant")])
        cols = [c for c in ["Species name", "Symbol", "substrate", "sub_CID",
                            "sub_SMILES", "Transformations", "PMID"] if c in pl.columns]
        pl[cols].drop_duplicates().to_csv(
            os.path.join(out, "plant_metabolite_inventory.tsv"), sep="\t", index=False)
        print(f"\n=== Plant metabolite inventory: {pl.substrate.nunique():,} compounds "
              f"from {pl['Species name'].nunique()} plants ===")
        print("    Join this to your host-plant table: it is the chemistry your")
        print("    herbivores actually face, with structures attached.")
    log(f"wrote {out}/")


def melt_reactions(rx):
    """One row per (enzyme, substrate). P450RDB stores up to four substrates
    per reaction in parallel column blocks."""
    keep = ["Symbol", "Name", "Uniprot ID", "EC number", "Species name",
            "Species", "Txid", "Transformations", "bond", "PMID", "DOI", "sequence"]
    keep = [c for c in keep if c in rx.columns]
    out = []
    for i in (1, 2, 3, 4):
        sc, cc = f"Substrate{i}", f"sub_CID{i}"
        sm = f"sub_Smiles{i}"
        if sc not in rx.columns:
            continue
        blk = rx[keep + [c for c in (sc, cc, sm) if c in rx.columns]].copy()
        blk = blk.rename(columns={sc: "substrate", cc: "sub_CID", sm: "sub_SMILES"})
        blk["substrate_slot"] = i
        out.append(blk)
    m = pd.concat(out, ignore_index=True)
    m = m[m.substrate.notna() & (m.substrate.astype(str).str.strip() != "")]
    # cofactors and co-substrates are not the chemistry you are modelling
    cof = {"o2", "h2o", "h+", "nadph", "nadp+", "nadh", "nad+", "fmn", "fmnh2",
           "fad", "fadh2", "co2", "glutathione", "electron"}
    m = m[~m.substrate.astype(str).str.strip().str.lower().isin(cof)]
    return m


def cmd_chem_triage(a):
    """Triage P450RDB before committing to it."""
    need_pandas()
    triage_p450rdb(a.dir, a.out)


# --------------------------------------------------------------------------- #
#  ICPD reference-table evidence triage
# --------------------------------------------------------------------------- #
# Evidence tiers, strongest first. A gene inherits the strongest tier of any
# paper that mentions it. The distinction is not cosmetic: overexpression in a
# resistant strain says nothing about which compound an enzyme turns over, and
# a substrate label built from correlative records will not survive review.
EVIDENCE_TIERS = [
    ("FUNCTIONAL_METABOLISM", [
        r"baculovirus", r"heterologously\s+express", r"recombinant\s+(protein|enzyme|P450|CYP)",
        r"metaboli[sz]ed?\s+(by|the|six|two|imperatorin|xanthotoxin)", r"metabolic\s+assay",
        r"pmol/min", r"\bin\s+vitro\s+metabolis", r"\bSf9\b", r"expressed\s+in\s+(E\.\s?coli|yeast|insect cells)",
        r"catalyz\w+\s+the", r"enzyme\s+converts", r"substrate\s+specificit",
        r"functional\s+characteri[sz]", r"enzyme\s+activit", r"kinetic\s+parameter",
        r"\bK\s?m\b", r"\bkcat\b", r"purified\s+(protein|enzyme)",
        r"microsome[s]?\s+(from|of)", r"turnover", r"depletion\s+assay",
        r"metaboli[sz]ing\s+activity", r"oxidat\w+\s+of\s+\w+"]),
    ("GENETIC_PERTURBATION", [
        r"RNA\s?i\b", r"RNA\s+interference", r"knock-?down", r"silenc", r"CRISPR",
        r"transgenic", r"piggyBac", r"mutant", r"loss[- ]of[- ]function",
        r"gene\s+editing", r"dsRNA", r"antisense"]),
    ("CORRELATIVE_EXPRESSION", [
        r"over-?express", r"up-?regulat", r"induc", r"expression\s+(profile|level|analysis)",
        r"qRT-?PCR", r"resistant\s+strain", r"microarray", r"transcriptom"]),
]
# Compounds worth pulling out of abstracts: these become the candidate ligand
# list for anything downstream.
# The functional insect-P450 literature is dominated by insecticide resistance
# in crop pests. In ICPD the insecticide:allelochemical ratio of gene-records is
# roughly 12:1. If the question is host-plant chemistry in non-pest herbivores,
# that is a second domain shift on top of any human-to-insect one, and it has to
# be stated rather than absorbed silently.
COMPOUND_CLASS = {
    "insecticide": ["imidacloprid", "thiamethoxam", "permethrin", "deltamethrin",
                    "cypermethrin", "cyhalothrin", "chlorantraniliprole", "phoxim",
                    "abamectin", "indoxacarb", "metaflumizone", "diazinon", "DDT",
                    "aldrin", "pyrethroid", "neonicotinoid", "organophosphate"],
    "synergist": ["piperonyl butoxide"],
    "allelochemical": ["imperatorin", "xanthotoxin", "furanocoumarin", "quercetin",
                       "flavone", "chlorogenic acid", "indole-3-carbinol", "rutin",
                       "tannin", "gossypol", "caffeine", "nicotine", "capsaicin",
                       "myristicin", "aflatoxin"],
    "endogenous": ["juvenile hormone", "ecdysone", "20-hydroxyecdysone",
                   "methyl farnesoate", "theobromine"],
}
COMPOUND_HINTS = [
    "imperatorin", "xanthotoxin", "furanocoumarin", "quercetin", "flavone",
    "chlorogenic acid", "indole-3-carbinol", "rutin", "tannin", "gossypol",
    "caffeine", "theobromine", "nicotine", "capsaicin", "myristicin",
    "piperonyl butoxide", "diazinon", "cypermethrin", "aldrin", "permethrin",
    "deltamethrin", "cyhalothrin", "phoxim", "imidacloprid", "thiamethoxam",
    "chlorantraniliprole", "abamectin", "indoxacarb", "metaflumizone",
    "DDT", "pyrethroid", "neonicotinoid", "organophosphate", "aflatoxin",
    "juvenile hormone", "ecdysone", "20-hydroxyecdysone", "methyl farnesoate",
]
GENE_SPLIT = re.compile(r"[;,]|\s{2,}")
CYP_NAME = re.compile(r"\bCYP\s?[0-9]+[A-Z]*[0-9]*(?:[-_v][0-9A-Za-z]+)?\b", re.I)


def cmd_icpd_evidence(a):
    """Sort ICPD reference records into evidence tiers.

    The References tables list a paper per row and often dozens of genes per
    paper. This explodes them to one row per (gene, paper), classifies each
    paper by the strongest evidence its abstract describes, and reports how
    many genes survive at each tier. Only FUNCTIONAL_METABOLISM records carry
    a defensible substrate label."""
    need_pandas()
    # Detect columns PER FILE. References1 names its gene column "Gene name or
    # symbol" and References2 names it "Gene symbol"; concatenating first and
    # detecting once picks whichever file came first and silently discards
    # every row of the other.
    tables = []
    for f in a.tables:
        if not os.path.exists(f):
            sys.exit(f"ERROR: not found: {f}")
        t = (pd.read_excel(f) if f.endswith((".xlsx", ".xls"))
             else pd.read_csv(f, sep=None, engine="python"))
        t["__source"] = os.path.basename(f)
        cols = {
            "gene": next((c for c in t.columns if re.search(r"gene", str(c), re.I)), None),
            "abstract": next((c for c in t.columns if re.search(r"abstract", str(c), re.I)), None),
            "title": next((c for c in t.columns if re.search(r"title", str(c), re.I)), None),
            "doi": next((c for c in t.columns if re.search(r"doi|reference", str(c), re.I)), None),
            "phenotype": next((c for c in t.columns if re.search(r"phenotype", str(c), re.I)), None),
            "effect": next((c for c in t.columns if re.search(r"variant.*effect|effect", str(c), re.I)), None),
            "species": next((c for c in t.columns if re.search(r"species", str(c), re.I)), None),
            "desc": next((c for c in t.columns if re.search(r"description", str(c), re.I)), None),
        }
        print(f"{os.path.basename(f)}: {t.shape[0]} rows")
        print("   " + "  ".join(f"{k}='{v}'" for k, v in cols.items() if v))
        if not cols["gene"]:
            print(f"   WARNING: no gene column; skipping this file")
            continue
        tables.append((t, cols))
    if not tables:
        sys.exit("ERROR: no usable table")

    rows = []
    for t, cols in tables:
      gcol, acol, tcol = cols["gene"], cols["abstract"], cols["title"]
      dcol, pcol = cols["doi"], cols["phenotype"]
      for r in t.itertuples(index=False):
        d = dict(zip(t.columns, r))
        text = " ".join(str(d.get(c, "")) for c in
                        (acol, tcol, cols["phenotype"], cols["effect"], cols["desc"]) if c)
        tier, hits = "UNCLASSIFIED", []
        for name, pats in EVIDENCE_TIERS:
            h = [p for p in pats if re.search(p, text, re.I)]
            if h:
                tier, hits = name, h
                break
        cpds = sorted({c for c in COMPOUND_HINTS if re.search(re.escape(c), text, re.I)})
        klass = sorted({k for k, v in COMPOUND_CLASS.items() if any(c in v for c in cpds)})
        genes = set()
        for chunk in GENE_SPLIT.split(str(d.get(gcol, ""))):
            genes.update(m.group(0).replace(" ", "").upper() for m in CYP_NAME.finditer(chunk))
            for w in chunk.split():
                w = w.strip("()[]").strip()
                if w and not CYP_NAME.match(w) and re.match(r"^[A-Za-z][A-Za-z0-9]{3,}$", w) \
                        and w.lower() not in ("none", "gene", "name", "symbol"):
                    genes.add(w)
        for g in sorted(genes):
            rows.append({"gene": g, "evidence_tier": tier,
                         "n_genes_in_paper": len(genes),
                         # A compound named in a paper covering 45 genes is
                         # associated with the PAPER, not with any one gene.
                         # Only single- or two-gene papers give a usable
                         # gene-compound pair without reading the paper.
                         "compound_confidence": ("HIGH" if len(genes) <= 2 else
                                                 "MEDIUM" if len(genes) <= 5 else "LOW"),
                         "species": str(d.get(cols["species"], "")) if cols["species"] else "",
                         "compounds_mentioned": "; ".join(cpds),
                         "compound_classes": "; ".join(klass),
                         "matched_patterns": "; ".join(hits[:3]),
                         "title": str(d.get(tcol, ""))[:200] if tcol else "",
                         "doi": str(d.get(dcol, "")) if dcol else "",
                         "phenotype": str(d.get(pcol, "")) if pcol else "",
                         "source_file": d.get("__source", "")})
      # end per-table loop
    ev = pd.DataFrame(rows)
    if ev.empty:
        sys.exit("ERROR: no gene names could be parsed")
    os.makedirs(a.out, exist_ok=True)
    ev.to_csv(os.path.join(a.out, "gene_paper_evidence.tsv"), sep="\t", index=False)

    order = {n: i for i, (n, _) in enumerate(EVIDENCE_TIERS)}
    order["UNCLASSIFIED"] = 99
    best = (ev.assign(_o=ev.evidence_tier.map(order))
              .sort_values("_o").groupby("gene").first().reset_index()
              .drop(columns=["_o"]))
    best.to_csv(os.path.join(a.out, "gene_best_evidence.tsv"), sep="\t", index=False)

    print("\n=== Papers by strongest evidence described ===")
    print(ev.drop_duplicates(["doi", "title"]).evidence_tier.value_counts().to_string())
    print("\n=== Distinct genes by strongest evidence anywhere ===")
    print(best.evidence_tier.value_counts().to_string())

    fn = best[best.evidence_tier == "FUNCTIONAL_METABOLISM"]
    print(f"\n=== {len(fn)} genes with a functional metabolism record ===")
    gold = fn[(fn.compound_confidence == "HIGH") & (fn.compounds_mentioned != "")]
    print(f"    of which {len(gold)} come from a paper covering at most two genes,")
    print(f"    so the gene-compound pairing can be trusted without reading it.")
    if len(gold):
        print("\n=== The set worth using as an external test set ===")
        print(gold[["gene", "species", "compounds_mentioned"]].head(30).to_string(index=False))
        gold.to_csv(os.path.join(a.out, "gene_compound_highconfidence.tsv"),
                    sep="\t", index=False)
    print("\nCompounds are matched at PAPER level. A paper on 45 genes naming five")
    print("compounds yields 225 gene-compound pairs and almost all are spurious:")
    print("that is why compound_confidence exists. Anything below HIGH needs the")
    print("paper read before it can carry a substrate label.")
    print("\nEverything else is correlative or genetic. Overexpression in a")
    print("resistant strain is not evidence that an enzyme metabolises a compound,")
    print("and RNAi knockdown is genetic evidence of involvement, not of substrate")
    print("specificity. Only the functional set can carry a substrate label.")

    big = ev.drop_duplicates(["doi", "title"]).nlargest(5, "n_genes_in_paper")
    if len(big):
        print("\n=== Papers listing the most genes (these are almost always expression screens) ===")
        print(big[["n_genes_in_paper", "evidence_tier", "title"]].to_string(index=False))

    kc = {}
    for ks in ev.compound_classes.dropna():
        for k in [x.strip() for x in ks.split(";") if x.strip()]:
            kc[k] = kc.get(k, 0) + 1
    if kc:
        print("\n=== What kind of chemistry has the literature actually tested? ===")
        for k, v in sorted(kc.items(), key=lambda x: -x[1]):
            print(f"  {k:16s} {v:6d} gene-records")
        if kc.get("insecticide", 0) > 3 * kc.get("allelochemical", 1):
            print("\n  The functional record is dominated by synthetic insecticides in")
            print("  crop pests. If your question is host-plant allelochemistry in")
            print("  non-pest herbivores, that is a domain shift you must report, not")
            print("  a set you can train on and quietly generalise from.")

    cnt = {}
    for cs in ev.compounds_mentioned.dropna():
        for c in [x.strip() for x in cs.split(";") if x.strip()]:
            cnt[c] = cnt.get(c, 0) + 1
    if cnt:
        cp = pd.DataFrame(sorted(cnt.items(), key=lambda x: -x[1]),
                          columns=["compound", "n_gene_records"])
        cp.to_csv(os.path.join(a.out, "compounds_mentioned.tsv"), sep="\t", index=False)
        print(f"\n=== Candidate compounds named across the abstracts: {len(cp)} ===")
        print(cp.head(15).to_string(index=False))
    log(f"wrote {a.out}/")


def cmd_cyp_benchmark(a):
    """Measure clan-assignment accuracy against ICPD's own curated clade calls.

    A large bitscore margin says an assignment is unambiguous given the four
    HMMs. It says nothing about whether the assignment is right. ICPD carries a
    curated Clade for every one of its 66,477 P450s, and the reference used only
    a few hundred per clan, so the remainder is a large held-out labelled test
    set that costs nothing to use.

    Sequences that went into the reference are excluded, so this is a genuine
    hold-out and not an accuracy estimate on the training data."""
    need_pandas()
    for exe in ("hmmsearch",):
        if not shutil.which(exe):
            sys.exit(f"ERROR: {exe} not found. module load HMMER")
    hmms = sorted(glob.glob(os.path.join(a.clan_hmms, "*.hmm")))
    if not hmms:
        sys.exit(f"ERROR: no .hmm files in {a.clan_hmms}\n"
                 "  build them first: defensome.py cyp --out <results> --refs <ref.faa>")
    os.makedirs(a.out, exist_ok=True)

    t = (pd.read_excel(a.table) if str(a.table).endswith((".xlsx", ".xls"))
         else pd.read_csv(a.table, sep=None, engine="python"))
    fasta_ids = {h.split()[0] for h, _ in read_fasta(a.fasta)}
    idc, best = None, 0
    for c in t.columns:
        n = len(fasta_ids & set(t[c].astype(str).str.strip()))
        if n > best:
            idc, best = c, n
    clc = next((c for c in t.columns if re.search(r"clan|clade", str(c), re.I)), None)
    ordc = next((c for c in t.columns if re.search(r"order", str(c), re.I)), None)
    cmp_ = next((c for c in t.columns if re.search(r"complete", str(c), re.I)), None)
    spc = next((c for c in t.columns if re.search(r"species", str(c), re.I)), None)
    if not (idc and clc):
        sys.exit(f"ERROR: need an ID and a clade column. Columns: {list(t.columns)}")
    log(f"truth from '{clc}', keyed on '{idc}' ({best}/{len(fasta_ids)} matched)")

    truth, meta = {}, {}
    keep_ord = {o.strip().lower() for o in (a.orders or "").split(",") if o.strip()}
    for r in t.itertuples(index=False):
        d = dict(zip(t.columns, r))
        sid = str(d[idc]).strip()
        cl = norm_clan(d.get(clc))
        if not cl:
            continue
        if keep_ord and ordc and not any(o in str(d.get(ordc, "")).lower() for o in keep_ord):
            continue
        if cmp_ and not str(d.get(cmp_, "")).strip().lower().startswith("full"):
            continue
        truth[sid] = cl
        meta[sid] = {"species": str(d.get(spc, "")) if spc else "",
                     "order": str(d.get(ordc, "")) if ordc else ""}

    used = set()
    if a.refs and os.path.exists(a.refs):
        used = {h.split()[0] for h, _ in read_fasta(a.refs)}
        log(f"excluding {len(used)} sequences that went into the reference")

    test = os.path.join(a.out, "heldout.faa")
    n_test, tested_ids, n_short, n_stop = 0, set(), 0, 0
    with open(test, "w") as fh:
        for hdr, seq in read_fasta(a.fasta):
            sid = hdr.split()[0]
            if sid not in truth or sid in used:
                continue
            if len(seq) < a.min_len:
                n_short += 1; continue
            if "*" in seq.rstrip("*"):
                n_stop += 1; continue
            fh.write(f">{sid}\n{seq.rstrip('*')}\n")
            tested_ids.add(sid); n_test += 1
    if not n_test:
        sys.exit("ERROR: held-out set is empty")
    log(f"held-out labelled sequences: {n_test} tested "
        f"({n_short} below --min-len {a.min_len}, {n_stop} with internal stops, excluded)")

    scores = {}
    for hmm in hmms:
        clan = os.path.basename(hmm)[:-4]
        tbl = os.path.join(a.out, f"hits_{clan}.tbl")
        run(["hmmsearch", "--noali", "-E", "1e-5", "--cpu", str(a.threads),
             "--tblout", tbl, "-o", os.devnull, hmm, test])
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                scores.setdefault(f[0], {})[clan] = float(f[5])

    rows = []
    for sid in truth:
        if sid in used or sid not in scores:
            continue
        sc = sorted(scores[sid].items(), key=lambda x: -x[1])
        pred, top = sc[0]
        margin = top - sc[1][1] if len(sc) > 1 else top
        rows.append((sid, meta[sid]["species"], truth[sid], pred,
                     round(top, 1), round(margin, 1), truth[sid] == pred))
    res = pd.DataFrame(rows, columns=["id", "species", "icpd_clade", "predicted",
                                      "score", "margin", "correct"])
    res.to_csv(os.path.join(a.out, "benchmark_calls.tsv"), sep="\t", index=False)
    unscored = len(tested_ids - set(scores))

    acc = res.correct.mean()
    print(f"\n=== Hold-out accuracy against ICPD's curated clade ===")
    print(f"  tested {len(res)} sequences; {unscored} of the tested set matched "
          f"no clan HMM at E<1e-5")
    print(f"  accuracy {acc*100:.2f}%  ({int(res.correct.sum())} correct, "
          f"{int((~res.correct).sum())} wrong)")
    cm = pd.crosstab(res.icpd_clade, res.predicted)
    print("\n=== Confusion matrix (rows = ICPD, columns = predicted) ===")
    print(cm.to_string())
    cm.to_csv(os.path.join(a.out, "confusion_matrix.tsv"), sep="\t")
    print("\n=== Per-clan recall and precision ===")
    for c in sorted(set(res.icpd_clade) | set(res.predicted)):
        tp = int(((res.icpd_clade == c) & (res.predicted == c)).sum())
        fn_ = int(((res.icpd_clade == c) & (res.predicted != c)).sum())
        fp = int(((res.icpd_clade != c) & (res.predicted == c)).sum())
        rec = tp / (tp + fn_) if tp + fn_ else float("nan")
        pre = tp / (tp + fp) if tp + fp else float("nan")
        print(f"  {c:6s} recall {rec*100:6.2f}%  precision {pre*100:6.2f}%  (n={tp+fn_})")
    print("\n=== Margin distribution, correct vs incorrect ===")
    for lbl, sub in (("correct", res[res.correct]), ("incorrect", res[~res.correct])):
        if len(sub):
            print(f"  {lbl:9s} n={len(sub):6d}  min {sub.margin.min():8.1f}  "
                  f"median {sub.margin.median():8.1f}  max {sub.margin.max():8.1f}")
    print("\n=== Accuracy against retention, by margin threshold ===")
    print("  threshold   kept        accuracy   errors kept")
    for th in (0, 10, 20, 50, 100, 117, 200, 300):
        k = res[res.margin >= th]
        if not len(k):
            continue
        print(f"  {th:9d}  {len(k):6d} ({len(k)/len(res)*100:5.1f}%)  "
              f"{k.correct.mean()*100:7.3f}%  {int((~k.correct).sum()):5d}")
    print("\n  Pick a threshold here rather than assuming one. If the minimum")
    print("  margin in your own dataset sits above the point where errors")
    print("  disappear, that is quantitative grounds for trusting every call.")

    print("\n=== Caveat that belongs in the results ===")
    print("  This compares two automated assignments. ICPD's clades come from")
    print("  its own pipeline, not from manual curation of every gene, so a")
    print("  disagreement identifies a gene where two methods differ, not a")
    print("  gene where yours is wrong. Report it as agreement, not accuracy.")

    worst = None
    try:
        per = [(c, ((res.icpd_clade == c) & (res.predicted == c)).sum() /
                max(int((res.icpd_clade == c).sum()), 1))
               for c in sorted(set(res.icpd_clade))]
        worst = min(per, key=lambda x: x[1])
    except Exception:
        pass
    if worst and worst[1] < 0.999:
        print(f"\n  Weakest clan: {worst[0]} at {worst[1]*100:.2f}% recall. If that is")
        print(f"  also your smallest reference, rerun cyp-refs with a larger")
        print(f"  --max-per-clan and rebuild; the two are usually connected.")

    wrong = res[~res.correct]
    if len(wrong):
        wrong.to_csv(os.path.join(a.out, "misassigned.tsv"), sep="\t", index=False)
        print(f"\n=== The {len(wrong)} disagreements ===")
        print(wrong.nlargest(min(15, len(wrong)), "margin")[
            ["id", "species", "icpd_clade", "predicted", "margin"]].to_string(index=False))
        print("\nHigh-margin disagreements are the interesting ones: either ICPD's")
        print("curated clade is wrong for that gene, or the clan boundary is not")
        print("where one of you thinks it is. Neither is a pipeline bug, and both")
        print("are worth a sentence in the results.")
    log(f"wrote {a.out}/")


# =========================================================================== #
#  Sample sheets, isoform collapsing, and method comparison
# =========================================================================== #
# A transcriptome-derived proteome is not comparable to a genome annotation
# until isoforms are collapsed. Trinity/EvidentialGene IDs carry an explicit
# isoform field (c0_g1_i21 is the 21st isoform of one gene), so counting
# proteins counts isoforms. Genome annotations vary: some headers carry a
# gene: tag, some use .t1 transcript suffixes, some are already gene-level.
# Getting this wrong inflates one arm of a comparison and invalidates it.
GENE_PATTERNS = [
    ("trinity_isoform", r"_i\d+\w*$",
     "Trinity/EvidentialGene isoform suffix (c0_g1_i1)"),
    ("transcript_t", r"\.t\d+$", "AUGUSTUS/BRAKER transcript suffix (.t1)"),
    ("transcript_p", r"\.p\d+$", "TransDecoder ORF suffix (.p1)"),
    ("transcript_dash_R", r"-R[A-Z]+\d*$", "FlyBase-style transcript suffix (-RA)"),
    ("version_dot", r"\.\d+$", "trailing version number, usually NOT an isoform"),
]
GENE_TAG = re.compile(r"\bgene[:=](\S+)")


def gene_of(header, regex=None):
    """Gene ID for a protein header. A gene: tag in the description always
    wins; otherwise the supplied regex is stripped from the first token."""
    m = GENE_TAG.search(header)
    if m:
        return m.group(1)
    tok = header.split()[0]
    return re.sub(regex, "", tok) if regex else tok


# "DToL" is not one annotation method. In practice a DToL release is BRAKER
# for some species, Ensembl for others, and occasionally raw AUGUSTUS. They
# differ systematically: on one ten-species set, BRAKER gave a median of about
# 18,900 genes and Ensembl about 14,000 after collapsing. Lumping them makes the
# reference arm of any comparison heterogeneous without anyone noticing.
_PIPELINE_HINTS = [("BRAKER", r"^BRAKER"), ("Ensembl", r"^ENS[A-Z]*P\d"),
                   ("AUGUSTUS", r"^AUGUSTUS"), ("Trinity", r"_c\d+_g\d+_i\d+"),
                   ("FlyBase", r"^FBpp"), ("RefSeq", r"^[XN]P_\d"),
                   ("TransDecoder", r"\.p\d+$")]


def _guess_pipeline(ids):
    for name, rx in _PIPELINE_HINTS:
        if sum(1 for i in ids if re.search(rx, i)) > len(ids) * 0.5:
            return name
    return "unknown"


_METHOD_HINTS = [("DToL", r"dtol|_pep_"),
                 ("genome_guided", r"_GG_|_GG$|genome.?guided"),
                 ("denovo", r"denovo|de.?novo|dnovo")]


def _guess_species_method(sample):
    meth = ""
    for name, rx in _METHOD_HINTS:
        if re.search(rx, sample, re.I):
            meth = name
            break
    if not meth:
        meth = "denovo"
    sp = re.sub(r"_GG(_\d{6,8})?", "", sample)
    sp = re.sub(r"(_\d{6,8})?(_okay.*|_pep.*|_DToL.*)$", "", sp, flags=re.I)
    sp = re.sub(r"_(GG|dtol|denovo|dnovo)$", "", sp, flags=re.I)
    return sp, meth


def cmd_inspect(a):
    """Examine every proteome, propose a gene-ID rule, write a draft sheet.

    Every failure in this pipeline has come from an ID assumption that was
    wrong for one input. This makes the assumption visible before anything is
    computed: it prints real headers, shows how many proteins collapse to how
    many genes under each candidate rule, and writes a sample sheet you edit
    rather than a guess you never see."""
    need_pandas()
    pats = ("*.fa", "*.faa", "*.fasta", "*.fa.gz", "*.faa.gz", "*.fasta.gz")
    files = sorted({f for p in pats for f in glob.glob(os.path.join(a.proteomes, p))})
    if not files:
        sys.exit(f"ERROR: no FASTA files in {a.proteomes}")
    rows = []
    for f in files:
        # Rule DETECTION uses a sample; the collapse ESTIMATE uses every header.
        # Sampling the first N for the estimate is biased: Trinity/EvidentialGene
        # output lists isoform-rich clusters first, which overstated transcriptome
        # collapse by up to 16 percentage points on real data.
        hdrs, allh, n, lens = [], [], 0, []
        for h, sq in read_fasta(f):
            n += 1
            allh.append(h)
            lens.append(len(sq.rstrip("*")))
            if len(hdrs) < 500:
                hdrs.append(h)
        base = os.path.basename(f)
        print(f"\n=== {base}  ({n:,} sequences) ===")
        for h in hdrs[:3]:
            print(f"    {h[:150]}")
        if sum(1 for h in hdrs if GENE_TAG.search(h)) > len(hdrs) * 0.8:
            rule, rx, desc = "gene_tag", "", "gene: tag present in the description"
            ng = len({gene_of(h) for h in hdrs})
        else:
            cands = []
            for name, crx, cdesc in GENE_PATTERNS:
                hit = sum(1 for h in hdrs if re.search(crx, h.split()[0]))
                if hit > len(hdrs) * 0.5:
                    ngx = len({gene_of(h, crx) for h in hdrs})
                    cands.append((name, crx, cdesc, hit, ngx))
                    print(f"    candidate {name:18s} matches {hit}/{len(hdrs)} -> "
                          f"{ngx} genes from {len(hdrs)} proteins")
            merging = [c for c in cands if c[4] < len(hdrs)]
            pick = (merging or cands)
            if pick:
                rule, rx, desc, _, ng = pick[0]
            else:
                rule, rx, desc = "none", "", "IDs already look gene-level"
                ng = len({h.split()[0] for h in hdrs})
        ng_all = len({gene_of(h, rx or None) for h in allh}) if rule != "none" \
            else len({h.split()[0] for h in allh})
        ratio = ng_all / n if n else 1.0
        print(f"    -> rule '{rule}': {desc}")
        print(f"    -> all {n:,} proteins collapse to {ng_all:,} genes "
              f"({(1-ratio)*100:.1f}%)")
        if ratio > 0.99 and re.search(r"_i\d", " ".join(hdrs[:50])):
            print("    WARNING: an isoform-looking suffix is present but nothing "
                  "collapsed. Check by hand.")
        utr = sum(1 for h in hdrs if "utrorf" in h.lower())
        if utr:
            print(f"    NOTE: {utr}/{len(hdrs)} headers carry 'utrorf'. EvidentialGene")
            print("          marks ORFs found in another transcript's UTR; many are")
            print("          spurious. Consider --drop-utrorf when collapsing.")
        sample = re.sub(r"\.gz$", "", base)
        sample = re.sub(r"\.(aa\.)?(fa|faa|fasta)$", "", sample)
        sp, meth = _guess_species_method(sample)
        pipe = _guess_pipeline([h.split()[0] for h in hdrs])
        est = ng_all
        # an insect genome annotation should land roughly 10k-30k genes
        flag = ""
        if meth == "DToL" and est > 35000:
            flag = "IMPLAUSIBLE_GENE_COUNT"
            print(f"    WARNING: {est:,} genes after collapsing is not plausible for an")
            print(f"             insect genome annotation (expect roughly 10k-30k).")
            if rule == "transcript_t":
                print("             The .t suffix only groups isoforms when the base is a")
                print("             GENE id (g123.t1, g123.t2). Here the base looks like a")
                print("             unique protein serial, so nothing groups. This is")
                print("             probably unfiltered ab initio output. See UPGRADING.md.")
        if pipe == "Trinity":
            print("    NOTE: Trinity 'genes' (c_g components) are not biological genes.")
            print("          A fragmented assembly splits one real gene across several")
            print("          components, so collapsed counts still run above a genome")
            print("          annotation. Use `compare --complete-only`.")
        shortest = min(lens) if lens else 0
        n_lt100 = sum(1 for x in lens if x < 100)
        if shortest >= 95 and n > 1000:
            print(f"    NOTE: the shortest protein is {shortest} aa and none are under 95.")
            print("          An ORF caller dropped short ORFs before this file was made")
            print("          (TransDecoder's default is -m 100). Metallothioneins (40-64 aa)")
            print("          cannot be in it. Add a `nucleotide` column pointing at the")
            print("          transcripts and run `short-orfs` to recover them.")
            flag = (flag + ";" if flag else "") + "NO_SHORT_PROTEINS"
        rows.append({"sample_id": sample, "species": sp, "method": meth,
                     "pipeline": pipe, "fasta": os.path.abspath(f),
                     "gene_regex": rx, "rule": rule, "n_proteins": n,
                     "est_genes": est, "collapse_pct": round((1 - ratio) * 100, 1),
                     "shortest_aa": shortest, "n_under_100aa": n_lt100,
                     "nucleotide": "", "flag": flag})
    sh = pd.DataFrame(rows)
    out = a.samplesheet or os.path.join(a.out, "samplesheet.tsv")
    d = os.path.dirname(os.path.abspath(out))
    if d:
        os.makedirs(d, exist_ok=True)
    sh.to_csv(out, sep="\t", index=False)
    print(f"\n=== draft sample sheet: {out} ===")
    print(sh[["sample_id", "species", "method", "pipeline", "rule", "n_proteins",
              "est_genes", "collapse_pct", "shortest_aa", "flag"]].to_string(index=False))
    if "pipeline" in sh and sh[sh.method == "DToL"].pipeline.nunique() > 1:
        print("\nNOTE: the DToL arm mixes annotation pipelines:")
        print(sh[sh.method == "DToL"].pipeline.value_counts().to_string())
        print("Recovery against 'DToL' will be recovery against different pipelines")
        print("for different species. compare reports it split by pipeline.")
    if sh.species.nunique() and sh.method.nunique() > 1:
        print("\n=== design ===")
        print(pd.crosstab(sh.species, sh.method).to_string())
    print("\nEDIT THIS FILE. species and method are guessed from the file name and")
    print("gene_regex from the headers. Everything downstream trusts it.")


def cmd_collapse(a):
    """Longest protein per gene, per the sample sheet, into one directory."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    for c in ("sample_id", "fasta"):
        if c not in sh.columns:
            sys.exit(f"ERROR: sample sheet needs a '{c}' column")
    missing = [r.fasta for r in sh.itertuples(index=False) if not os.path.exists(r.fasta)]
    if missing:
        sys.exit("ERROR: missing FASTA files:\n  " + "\n  ".join(missing[:5]))
    dest = os.path.join(a.out, "proteomes")
    os.makedirs(dest, exist_ok=True)
    rows = []
    for r in sh.itertuples(index=False):
        rx = getattr(r, "gene_regex", "")
        rx = "" if rx is None or (isinstance(rx, float) and np.isnan(rx)) else str(rx).strip()
        best, n_in, n_utr = {}, 0, 0
        for h, seq in read_fasta(r.fasta):
            n_in += 1
            if a.drop_utrorf and "utrorf" in h.lower():
                n_utr += 1
                continue
            g = gene_of(h, rx or None)
            if g not in best or len(seq) > len(best[g][1]):
                best[g] = (h.split()[0], seq)
        with open(os.path.join(dest, f"{r.sample_id}.faa"), "w") as fh:
            for g, (pid, seq) in best.items():
                fh.write(f">{pid}\n{seq.rstrip('*')}\n")
        rows.append({"sample_id": r.sample_id,
                     "species": getattr(r, "species", ""),
                     "method": getattr(r, "method", ""),
                     "n_proteins_in": n_in, "n_genes_out": len(best),
                     "dropped_utrorf": n_utr,
                     "collapse_pct": round((1 - len(best) / max(n_in, 1)) * 100, 1)})
        log(f"{r.sample_id}: {n_in} proteins -> {len(best)} genes"
            + (f", {n_utr} utrorf dropped" if n_utr else ""))
        if getattr(r, "method", "") == "DToL" and len(best) > 35000:
            print(f"    WARNING: {r.sample_id} has {len(best):,} genes after collapsing;")
            print(f"             an insect genome should be roughly 10k-30k. Probably")
            print(f"             unfiltered ab initio output. See UPGRADING.md, or")
            print(f"             leave it out: compare --exclude {getattr(r, 'species', r.sample_id)}")
    st = pd.DataFrame(rows)
    st.to_csv(os.path.join(a.out, "collapse_stats.tsv"), sep="\t", index=False)
    print("\n" + st.to_string(index=False))
    if st.method.nunique() > 1:
        print("\n=== median by method ===")
        print(st.groupby("method")[["n_proteins_in", "n_genes_out", "collapse_pct"]]
                .median().round(1).to_string())
        print("\nIf one method collapses far more than another, its protein count was")
        print("measuring isoforms. Compare n_genes_out, never n_proteins_in.")
    log(f"collapsed proteomes in {dest}")


def cmd_compare(a):
    """Compare defensome recovery across annotation methods, paired by species.

    The design is what makes this work: the same species assayed by several
    methods turns annotation source from a confound into a measured variable,
    and pairs every comparison within species. Cross-sectional differences
    between 30 proteomes would tell you almost nothing; ten paired triplets
    tell you how much of a genome's defensome a transcriptome recovers."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    need = {"sample_id", "species", "method"}
    if not need.issubset(sh.columns):
        sys.exit(f"ERROR: sample sheet needs {sorted(need)}")
    cf = os.path.join(a.out, "counts.tsv")
    if not os.path.exists(cf):
        sys.exit(f"ERROR: {cf} not found. `scan` and `annotate` have not completed.\n"
                 "  python3 defensome.py doctor   checks the tools they need")
    counts = pd.read_csv(cf, sep="\t", index_col=0)
    npath = os.path.join(a.out, "counts_per10k.tsv")
    norm = pd.read_csv(npath, sep="\t", index_col=0) if os.path.exists(npath) else None
    dm = read_map(a.map)
    outdir = os.path.join(a.out, "comparison"); os.makedirs(outdir, exist_ok=True)

    if a.exclude:
        drop = {x.strip() for x in a.exclude.split(",") if x.strip()}
        log(f"excluding samples: {sorted(drop)}")
        sh = sh[~sh.sample_id.isin(drop) & ~sh.species.isin(drop)]
    sh = sh[sh.sample_id.isin(counts.index)]
    if sh.empty:
        sys.exit("ERROR: no sample_id in the sheet matches counts.tsv.\n"
                 "  The scan/annotate steps must use the collapsed proteomes, whose\n"
                 "  file names are the sample_id values.")
    meta = sh.set_index("sample_id")
    methods = sorted(meta.method.unique())
    species = sorted(meta.species.unique())
    design = pd.crosstab(meta.species, meta.method)
    print("=== design ===")
    print(design.to_string())
    complete = design[(design > 0).all(axis=1)].index.tolist()
    print(f"\n{len(complete)} of {len(species)} species have every method; "
          f"paired tests use those {len(complete)}.")
    if len(complete) < 3:
        print("WARNING: too few complete sets for a paired test; reporting "
              "descriptives only.")

    core = [f for f in dm.loc[dm.tier == "CORE", "family"]
            if f in counts.columns and counts[f].sum() > 0]
    # Within one species the true defensome is fixed, so for a PAIRED method
    # comparison raw counts are the right measure. Dividing by proteome size
    # re-injects the artefact: the denominator is annotation-dependent (Trinity
    # fragmentation inflates it, an unfiltered ab initio set inflates it more),
    # so per-10k would penalise exactly the arms with the noisiest gene totals.
    # Per-10k is for comparing DIFFERENT species; it is opt-in here.
    if a.complete_only:
        cf = os.path.join(a.out, "domains", "completeness.tsv")
        if not os.path.exists(cf):
            sys.exit(f"ERROR: --complete-only needs {cf}. Run `domains` first.")
        comp = pd.read_csv(cf, sep="\t")
        counts = (comp[comp.status == "COMPLETE"]
                  .groupby(["species", "family"]).protein.nunique()
                  .unstack(fill_value=0).reindex(index=counts.index,
                                                 columns=counts.columns, fill_value=0))
        log("counting COMPLETE domain architectures only (fragments excluded)")
    use = norm if (norm is not None and a.per10k and not a.complete_only) else counts
    lab = ("copies per 10k genes" if use is norm else
           "COMPLETE gene counts" if a.complete_only else "raw gene counts")

    # long table: one row per species x method x family
    rows = []
    for sid in counts.index:
        if sid not in meta.index:
            continue
        for f in core:
            rows.append({"species": meta.loc[sid, "species"],
                         "method": meta.loc[sid, "method"], "sample_id": sid,
                         "family": f, "raw": int(counts.loc[sid, f]),
                         "value": float(use.loc[sid, f])})
    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(outdir, "long_counts.tsv"), sep="\t", index=False)

    piv = long.pivot_table(index=["species", "family"], columns="method",
                           values="value").reset_index()
    piv.to_csv(os.path.join(outdir, "species_family_by_method.tsv"), sep="\t", index=False)

    print(f"\n=== Median {lab} by method, summed over CORE families ===")
    tot = long.groupby(["species", "method"]).value.sum().reset_index()
    print(tot.pivot(index="species", columns="method", values="value")
             .round(1).to_string())

    # recovery against the reference method
    ref = a.reference if a.reference in methods else (
        "DToL" if "DToL" in methods else methods[0])
    print(f"\n=== Recovery relative to '{ref}' (raw gene counts, paired within species) ===")
    rec_rows = []
    for sp_ in complete:
        base = long[(long.species == sp_) & (long.method == ref)].set_index("family").raw
        for m in methods:
            if m == ref:
                continue
            other = long[(long.species == sp_) & (long.method == m)].set_index("family").raw
            for f in core:
                b, o = base.get(f, 0), other.get(f, 0)
                rec_rows.append({"species": sp_, "method": m, "family": f,
                                 "ref_count": b, "count": o,
                                 "recovery": (o / b) if b else np.nan})
    rec = pd.DataFrame(rec_rows)
    if not rec.empty and "pipeline" in meta.columns:
        refpipe = (meta[meta.method == ref].reset_index()
                   .set_index("species").pipeline.to_dict())
        rec["ref_pipeline"] = rec.species.map(refpipe)
    if not rec.empty:
        rec.to_csv(os.path.join(outdir, "recovery.tsv"), sep="\t", index=False)
        summ = (rec.groupby("method")
                   .agg(median_recovery=("recovery", "median"),
                        families_lost=("recovery", lambda x: int((x == 0).sum())),
                        families_exceeding=("recovery", lambda x: int((x > 1).sum())),
                        n=("recovery", "size")).round(3))
        print(summ.to_string())
        if "ref_pipeline" in rec and rec.ref_pipeline.nunique() > 1:
            print(f"\n=== Same, split by the pipeline that produced the '{ref}' reference ===")
            print(rec.groupby(["ref_pipeline", "method"]).recovery.median()
                     .unstack().round(3).to_string())
            print("If these differ, the reference arm is not one thing and a single")
            print("'recovery against DToL' number averages over two denominators.")
        print("\nRecovery above 1 does not mean the transcriptome found extra genes.")
        print("It usually means residual isoform inflation, or a family where the")
        print("genome annotation missed members. Check those families by hand.")
        byfam = (rec.groupby(["method", "family"]).recovery.median()
                    .unstack(0).round(2).sort_values(rec.method.unique()[0]))
        byfam.to_csv(os.path.join(outdir, "recovery_by_family.tsv"), sep="\t")
        print(f"\n=== Worst-recovered CORE families ===")
        print(byfam.head(10).to_string())

    # paired tests
    if len(complete) >= 3:
        try:
            from scipy import stats as _st
        except ImportError:
            _st = None
        if _st:
            trows = []
            for f in core:
                for i, m1 in enumerate(methods):
                    for m2 in methods[i + 1:]:
                        x, y = [], []
                        for sp_ in complete:
                            v1 = long[(long.species == sp_) & (long.method == m1) &
                                      (long.family == f)].value
                            v2 = long[(long.species == sp_) & (long.method == m2) &
                                      (long.family == f)].value
                            if len(v1) and len(v2):
                                x.append(float(v1.iloc[0])); y.append(float(v2.iloc[0]))
                        if len(x) >= 3 and any(p != q for p, q in zip(x, y)):
                            try:
                                w = _st.wilcoxon(x, y)
                                trows.append({"family": f, "method_a": m1, "method_b": m2,
                                              "n_pairs": len(x),
                                              "median_a": float(np.median(x)),
                                              "median_b": float(np.median(y)),
                                              "W": float(w.statistic), "p": float(w.pvalue)})
                            except Exception:
                                pass
            if trows:
                tt = pd.DataFrame(trows).sort_values("p")
                tt["p_bh"] = tt.p * len(tt) / tt.p.rank()
                tt["p_bh"] = tt.p_bh[::-1].cummin()[::-1].clip(upper=1)
                tt.to_csv(os.path.join(outdir, "paired_tests.tsv"), sep="\t", index=False)
                print(f"\n=== Wilcoxon signed-rank, paired within species ({len(complete)} pairs) ===")
                print(tt.head(12).round(4).to_string(index=False))
                print(f"\n{int((tt.p_bh < 0.05).sum())} of {len(tt)} family-by-method-pair")
                print("comparisons survive a 5% false discovery rate. These are paired")
                print("within species, so unlike a cross-species diet test they are not")
                print("confounded by phylogeny: the same genome is being measured twice.")

    _compare_figures(outdir, long, rec, core, methods, complete, lab, ref)
    log(f"wrote {outdir}/")


def _compare_figures(outdir, long, rec, core, methods, complete, lab, ref):
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
    except ImportError:
        log("matplotlib unavailable, skipping comparison figures"); return
    _pub_style()
    fd = os.path.join(outdir, "figures")
    if os.path.exists(os.path.join(fd, "FIGURES.md")):
        os.remove(os.path.join(fd, "FIGURES.md"))
    W = FigureWriter(fd, "comparison_all_figures.pdf")
    methods = method_order(methods)
    cols = {m: method_colour(m, i) for i, m in enumerate(methods)}
    xs = {m: i for i, m in enumerate(methods)}

    # C1: slopegraph, one line per species across methods, y from zero
    tot = long.groupby(["species", "method"]).value.sum().unstack().reindex(columns=methods)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
    ends = []
    for sp_ in tot.index:
        v = tot.loc[sp_]
        pts = [(xs[m], v[m]) for m in methods if not np.isnan(v[m])]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#adb5bd", lw=1, zorder=1)
        for m in methods:
            if not np.isnan(v[m]):
                ax.scatter(xs[m], v[m], s=36, color=cols[m], zorder=3, edgecolor="white", lw=.6)
        last = [m for m in methods if not np.isnan(v[m])][-1]
        ends.append([v[last], xs[last], sp_])
    # push end labels apart so species with equal totals stay readable
    ends.sort()
    span = (np.nanmax(tot.values) or 1) * .045
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < span:
            ends[k][0] = ends[k - 1][0] + span
    for ylab, xpos, sp_ in ends:
        yreal = tot.loc[sp_, methods[int(xpos)]]
        ax.annotate(sp_.replace("_", " "), (xpos, yreal), xytext=(xpos + .08, ylab),
                    textcoords="data", va="center", fontsize=7.5, style="italic",
                    arrowprops=dict(arrowstyle="-", color="#ced4da", lw=.5)
                    if abs(ylab - yreal) > 1e-9 else None)
    ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods)
    ax.set_xlim(-.3, len(methods) - .3); ax.set_ylim(0, None)
    ax.set_ylabel(f"total CORE defensome ({lab})")
    ax.set_title("A  Same genome, three annotation routes", loc="left")
    if ref in tot.columns:
        others = [m for m in methods if m != ref]
        ratio = tot[others].div(tot[ref], axis=0)
        for i, m in enumerate(others):
            vals = ratio[m].dropna()
            ax2.scatter(np.random.default_rng(1).normal(i, .05, len(vals)), vals, s=30,
                        color=cols[m], edgecolor="white", lw=.6, zorder=3)
            if len(vals):
                ax2.hlines(vals.median(), i - .25, i + .25, color=cols[m], lw=2.5, zorder=4)
                ax2.annotate(f"median {vals.median():.2f}", (i + .28, vals.median()),
                             va="center", fontsize=8, color=cols[m])
        ax2.axhline(1, ls="--", color="#6c757d", lw=1)
        ax2.set_xticks(range(len(others))); ax2.set_xticklabels(others)
        ax2.set_xlim(-.5, len(others) - .1)
        ax2.set_ylabel(f"total relative to {ref}")
        ax2.set_title(f"B  Recovery of the {ref} total", loc="left")
    fig.tight_layout()
    W.save(fig, "C1_paired_totals",
           f"(A) Total CORE defensome per species under each annotation route, joined within "
           f"species; y-axis from zero. (B) Each route's total divided by the {ref} total for the "
           f"same species; dashed line is parity, bar is the median. Values are {lab}.")

    # C2: Cleveland dot plot, scales to any number of families
    med = long.groupby(["family", "method"]).value.median().unstack().reindex(columns=methods)
    sort_by = ref if ref in med.columns else methods[0]
    med = med.sort_values(sort_by, ascending=True)
    fig, ax = plt.subplots(figsize=(7.5, max(3.2, .32 * len(med) + 1.2)))
    dodge = {m: (i - (len(methods) - 1) / 2) * .18 for i, m in enumerate(methods)}
    for yi, (fam, r) in enumerate(med.iterrows()):
        v = r.dropna()
        if len(v):
            ax.plot([v.min(), v.max()], [yi, yi], color="#dee2e6", lw=3, zorder=1,
                    solid_capstyle="round")
        for m in methods:
            if not np.isnan(r[m]):
                # vertical dodge: equal values would otherwise hide all but one colour
                ax.scatter(r[m], yi + dodge[m], s=38, color=cols[m], zorder=3,
                           edgecolor="white", lw=.6, label=m if yi == 0 else None)
    ax.set_yticks(range(len(med))); ax.set_yticklabels(med.index)
    pos = med.values[~np.isnan(med.values)]
    if len(pos) and pos.max() > 0 and pos[pos > 0].min() > 0 and pos.max() / pos[pos > 0].min() > 25:
        ax.set_xscale("log")
    ax.set_xlabel(f"median across species ({lab})")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", ncol=len(methods))
    ax.set_title("Per-family gene counts by annotation route", loc="left")
    fig.tight_layout()
    W.save(fig, "C2_family_by_method",
           f"Median {lab} per family across species, one dot per annotation route. The grey bar "
           f"spans the routes, so a long bar is a family whose count depends heavily on how the "
           f"proteome was made. Families sorted by the {sort_by} value.")

    # C3: recovery heatmap, family x method, centred on parity
    if rec is not None and not rec.empty:
        others = [m for m in methods if m != ref]
        hm = rec.groupby(["family", "method"]).recovery.median().unstack().reindex(columns=others)
        hm = hm.loc[hm.mean(axis=1).sort_values().index]
        fig, ax = plt.subplots(figsize=(2.2 + 1.5 * len(others), max(3.2, .3 * len(hm) + 1.4)))
        data = hm.values.astype(float)
        vmax = max(2.0, np.nanmax(data) if np.isfinite(np.nanmax(data)) else 2.0)
        im = ax.imshow(np.nan_to_num(data, nan=1.0), cmap="RdBu", aspect="auto",
                       norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax))
        for yi in range(data.shape[0]):
            for xi in range(data.shape[1]):
                v = data[yi, xi]
                if np.isnan(v):
                    ax.add_patch(plt.Rectangle((xi - .5, yi - .5), 1, 1, color="#f1f3f5"))
                    ax.text(xi, yi, "n/a", ha="center", va="center", fontsize=7, color="#adb5bd")
                else:
                    ax.text(xi, yi, f"{v:.2f}", ha="center", va="center", fontsize=7.5,
                            color="white" if abs(v - 1) > .45 else "#212529")
        ax.set_xticks(range(len(others))); ax.set_xticklabels(others)
        ax.set_yticks(range(len(hm))); ax.set_yticklabels(hm.index)
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, shrink=.6, pad=.02)
        cb.set_label(f"median recovery vs {ref}")
        ax.set_title(f"Recovery relative to {ref}", loc="left")
        fig.tight_layout()
        W.save(fig, "C3_recovery_heatmap",
               f"Median over species of (gene count in route / gene count in {ref}) for each "
               f"family. 1 is parity; red below, blue above. Above 1 usually means residual "
               f"isoform or fragment inflation rather than genes {ref} missed. n/a: {ref} had none.")
    W.close()
    log(f"comparison figures in {fd} (PDF, PNG, combined PDF, FIGURES.md)")


def cmd_doctor(a):
    """Report which tools and libraries are present, and which steps can run."""
    import platform
    print(f"defensome {__version__}  ({os.path.abspath(__file__)})")
    print(f"python    {platform.python_version()}  ({sys.executable})\n")
    libs = {}
    for m in ("pandas", "numpy", "scipy", "matplotlib", "openpyxl"):
        try:
            mod = __import__(m)
            libs[m] = getattr(mod, "__version__", "yes")
        except Exception:
            libs[m] = None
    print("python libraries")
    for m, v in libs.items():
        need = "required" if m in ("pandas", "numpy") else "optional"
        print(f"  {'OK ' if v else '-- '} {m:11s} {v or 'missing':12s} {need}")
    print("\nexternal tools")
    tools = {}
    for t in ("hmmsearch", "hmmbuild", "mafft", "FastTree", "cd-hit", "node"):
        p = tool(t) or (tool("fasttree") if t == "FastTree" else None)
        ver = ""
        if p:
            try:
                o = subprocess.run([p, "-h"] if t != "node" else [p, "--version"],
                                   capture_output=True, text=True, timeout=10)
                txt = (o.stdout + o.stderr)
                m_ = re.search(r"(\d+\.\d+(?:\.\d+)?)", txt)
                ver = m_.group(1) if m_ else ""
            except Exception:
                pass
        tools[t] = p
        print(f"  {'OK ' if p else '-- '} {t:11s} {ver or ('missing' if not p else ''):12s} {p or ''}")
    steps = {
        "inspect / collapse / compare / report / dashboard": libs["pandas"] and libs["numpy"],
        "scan": tools["hmmsearch"],
        "cyp": tools["hmmbuild"] and tools["hmmsearch"] and tools["mafft"],
        "trees": tools["mafft"] and tools["FastTree"],
        "figures": libs["matplotlib"],
        "dashboard test": tools["node"],
    }
    print("\nsteps that can run here")
    for k, ok in steps.items():
        print(f"  {'OK ' if ok else '-- '} {k}")
    missing = [k for k, ok in steps.items() if not ok]
    if missing:
        print("\nOn an Lmod cluster, typically:")
        print("  module load HMMER MAFFT FastTree")
        print("If HMMER and your Python need incompatible toolchains, point at the")
        print("binary instead of loading both modules:")
        print("  module load HMMER && export DEFENSOME_HMMSEARCH=$(which hmmsearch) \\")
        print("                    && export DEFENSOME_HMMBUILD=$(which hmmbuild)")
        print("  module purge && module load Python")
    if a.require:
        req = {x.strip() for x in a.require.split(",") if x.strip()}
        bad = [r for r in req if r in ("scan",) and not tools["hmmsearch"]] + \
              [r for r in req if r in ("cyp",) and not steps["cyp"]] + \
              [r for r in req if r in ("trees",) and not steps["trees"]] + \
              [r for r in req if r in ("core",) and not steps[
                  "inspect / collapse / compare / report / dashboard"]]
        if bad:
            sys.exit(f"\nERROR: required step(s) cannot run: {', '.join(sorted(set(bad)))}")


# --------------------------------------------------------------------------- #
#  Short-ORF recovery from nucleotide sequence
# --------------------------------------------------------------------------- #
_COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")
_CODON = {}
for _i, _a in enumerate("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _CODON["TCAG"[_i // 16] + "TCAG"[(_i // 4) % 4] + "TCAG"[_i % 4]] = _a


def _translate(nt):
    return "".join(_CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def short_orfs(seq, min_aa=25, max_aa=99):
    """ATG-to-stop ORFs of min_aa..max_aa residues on both strands.

    Stop and ATG positions are found once with regular expressions and paired
    by frame, so only the candidate ORFs are ever translated. Translating whole
    transcriptomes six ways in pure Python would be far slower."""
    import bisect
    seq = seq.upper()
    out = []
    for strand, s in (("+", seq), ("-", seq.translate(_COMP)[::-1])):
        stops = [[], [], []]
        for m in re.finditer(r"(?=(TAA|TAG|TGA))", s):
            stops[m.start() % 3].append(m.start())
        atgs = [[], [], []]
        for m in re.finditer(r"(?=ATG)", s):
            atgs[m.start() % 3].append(m.start())
        for f in range(3):
            prev = f - 3
            for st in stops[f] + [None]:
                end = st if st is not None else len(s) - ((len(s) - f) % 3)
                j = bisect.bisect_left(atgs[f], prev + 3)
                if j < len(atgs[f]) and atgs[f][j] < end and st is not None:
                    a0 = atgs[f][j]
                    n_aa = (st - a0) // 3
                    if min_aa <= n_aa <= max_aa:
                        out.append((strand, a0, n_aa, _translate(s[a0:st])))
                if st is None:
                    break
                prev = st
    return out


def cmd_short_orfs(a):
    """Recover proteins too short for the ORF caller that built the proteome.

    TransDecoder keeps ORFs of at least 100 residues by default, and most gene
    predictors skip very short genes, so an insect metallothionein (40-64
    residues) is usually absent from a transcriptome proteome before any
    search runs. No domain rule can find a protein that is not in the file.
    This goes back to the transcripts named in the sample sheet's `nucleotide`
    column, takes ATG-to-stop ORFs in the band the ORF caller discarded, and
    applies the metallothionein composition rule. Genes, not transcripts, are
    counted, using the sample sheet's gene_regex."""
    need_pandas()
    sh = pd.read_csv(a.samplesheet, sep="\t")
    col = a.nucleotide_col
    if col not in sh.columns:
        sys.exit(f"ERROR: the sample sheet has no '{col}' column.\n"
                 f"  Add one giving the transcript FASTA for each sample (Trinity.fasta etc.).\n"
                 f"  Leave it empty for samples without one; they are skipped.")
    d = os.path.join(a.out, "short_orfs"); os.makedirs(d, exist_ok=True)
    hits, summ = [], []
    for r in sh.itertuples(index=False):
        nt = getattr(r, col, None)
        if nt is None or (isinstance(nt, float) and np.isnan(nt)) or not str(nt).strip():
            continue
        nt = str(nt).strip()
        if not os.path.exists(nt):
            log(f"WARNING: {r.sample_id}: nucleotide file not found: {nt}"); continue
        if os.path.getsize(nt) > a.max_mb * 1e6:
            log(f"WARNING: {r.sample_id}: {nt} is over {a.max_mb} MB; this step is meant for "
                f"transcript sets, not whole genomes. Use --max-mb to override.")
            continue
        rx = getattr(r, "gene_regex", "")
        rx = "" if rx is None or (isinstance(rx, float) and np.isnan(rx)) else str(rx).strip()
        n_tx = n_orf = 0
        genes = {}
        with open(os.path.join(d, f"{r.sample_id}.mt_like.faa"), "w") as fo:
            for hdr, seq in read_fasta(nt):
                n_tx += 1
                for strand, pos, n_aa, prot in short_orfs(seq, a.min_aa, a.max_aa):
                    n_orf += 1
                    if not _is_mt_like(prot):
                        continue
                    tid = hdr.split()[0]
                    gene = gene_of(hdr, rx or r"_i\d+\w*$")
                    f_ = mt_features(prot)
                    fo.write(f">{tid}|{strand}{pos} gene={gene} len={n_aa}\n{prot}\n")
                    hits.append({"sample_id": r.sample_id, "species": getattr(r, "species", ""),
                                 "method": getattr(r, "method", ""), "gene": gene,
                                 "transcript": tid, "strand": strand, "start": pos,
                                 **f_, "protein": prot})
                    genes.setdefault(gene, prot)
        summ.append({"sample_id": r.sample_id, "species": getattr(r, "species", ""),
                     "method": getattr(r, "method", ""), "transcripts": n_tx,
                     "short_orfs_scanned": n_orf, "mt_like_genes": len(genes)})
        log(f"{r.sample_id}: {n_tx:,} transcripts, {n_orf:,} short ORFs, "
            f"{len(genes)} MT-like genes")
    if not summ:
        sys.exit(f"ERROR: no sample had a usable '{col}' entry")
    hd = pd.DataFrame(hits)
    hd.to_csv(os.path.join(a.out, "short_orf_mt.tsv"), sep="\t", index=False)
    sd = pd.DataFrame(summ)
    sd.to_csv(os.path.join(a.out, "short_orf_summary.tsv"), sep="\t", index=False)
    print("\n" + sd.to_string(index=False))
    print("\nThese are composition-screen candidates from ORFs the proteome never")
    print("contained. They are reported apart from every count table. Confirm a few")
    print("by alignment to a known insect metallothionein before citing numbers.")
    log(f"wrote {a.out}/short_orf_mt.tsv and short_orf_summary.tsv")


def cmd_dataset_help(a):
    print(DATASET_HELP)


DATASET_HELP = """
External datasets: what each can and cannot do
==============================================

                what it holds                    protein?  negatives?  insect?
  ICPD          66,477 P450s / 682 species       yes       no          yes
  P450RDB       enzyme -> substrate + sequence   yes       NO          check
  Ni 2025       ~2000 compounds per human CYP    none      YES         no

1. ICPD  http://www.insectp450.net/ui/#/page/download

   Once downloaded (all_prot.fa, information.xlsx, References1.xlsx,
   References2.xlsx):

       python3 defensome.py cyp-refs --source icpd \\
           --fasta db/icpd/all_prot.fa --table db/icpd/information.xlsx \\
           --orders Lepidoptera --out db/icpd/ref.faa
       cd-hit -i db/icpd/ref.faa -o db/icpd/ref.c80.faa -c 0.8 -n 5
       python3 defensome.py cyp --out results/ --refs db/icpd/ref.c80.faa

       python3 defensome.py icpd-evidence \\
           --tables db/icpd/References1.xlsx db/icpd/References2.xlsx \\
           --out db/icpd/evidence

   The FASTA is keyed on the table's "Database ID" column (Abter001), NOT
   "Protein ID"; cyp-refs finds the right column by matching against the FASTA
   headers rather than guessing from the column name. Entries marked
   "fragmented" in the Completeness column are excluded by default: a truncated
   P450 degrades any profile it goes into.

   OLD INSTRUCTIONS BELOW, superseded:
   Download 'Protein Sequence' (fasta) and 'Sequence information table' (xlsx).
   The links are behind a web UI, so fetch them in a browser and copy across:

       mkdir -p db/icpd
       # then scp/rsync the two files into db/icpd/

   Build a lepidopteran clan reference and use it:

       python3 defensome.py cyp-refs --source icpd \
           --fasta db/icpd/ICPD_protein.fasta \
           --table db/icpd/ICPD_sequence_info.xlsx \
           --orders Lepidoptera --out db/icpd/ref.faa
       cd-hit -i db/icpd/ref.faa -o db/icpd/ref.c80.faa -c 0.8 -n 5
       python3 defensome.py cyp --out results/ --refs db/icpd/ref.c80.faa

   Validate the clan assignment against ICPD's own curated clades. The
   reference uses only a few hundred sequences per clan, so the rest is a large
   held-out labelled test set:

       python3 defensome.py cyp-benchmark \\
           --fasta db/icpd/all_prot.fa --table db/icpd/information.xlsx \\
           --clan-hmms results/cyp/clan_hmms --refs db/icpd/ref.c80.faa \\
           --orders Lepidoptera --out db/icpd/benchmark

   A large bitscore margin says an assignment is unambiguous. It does not say
   it is right. This measures accuracy directly and writes a confusion matrix.

   Also download 'References2' (genes associated with phenotypes). Separating
   'overexpressed in a resistant strain' (correlative) from 'heterologously
   expressed and shown to metabolise X' (functional) is the single most
   important curation decision in this analysis. Only the second is a
   substrate label.

   Cite: Wu et al. (2025) Mol Ecol Resour 25:e14070.

2. P450RDB  https://www.cellknowledge.com.cn/p450rdb_v2/download.html
       python3 defensome.py chem-triage --dir db/p450rdb --out db/p450rdb/triage

   Expect mostly plant and microbial BIOSYNTHETIC reactions. A reaction
   database has no negatives, so it cannot train a binary classifier. Its real
   value is triage/plant_metabolite_inventory.tsv: the chemistry your
   herbivores face, with structures, joinable to your host-plant table.

3. Ni et al. (2025) Sci Data 12:1427, doi 10.1038/s41597-025-05753-8
   The only source with curated NON-substrates. Two jobs, neither of them
   predicting insect biology: reproduce their GCN as a methods benchmark
   (published MCC 0.51 to 0.72), then measure how far a human-trained model
   transfers to functionally validated insect pairs. That drop is the result.
   Licence is CC BY-NC-ND; NoDerivatives constrains merged datasets. Check the
   Figshare record separately.

Nothing you already have under db/ can substitute for ICPD. Until it is
downloaded, --source flybase keeps you running with Drosophila references and
LOW-confidence calls.
"""


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(prog="defensome.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version",
                   version=f"defensome {__version__}  ({os.path.abspath(__file__)})")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(q, need_proteomes=False):
        q.add_argument("--out", required=True, help="output directory")
        q.add_argument("--map", default=DEFAULT_MAP, help="defensome map TSV")
        q.add_argument("--proteomes", required=need_proteomes,
                       help="directory of .faa proteomes")

    s = sub.add_parser("scan", help="run hmmsearch per proteome")
    common(s, True)
    s.add_argument("--pfam", required=True, help="Pfam-A.hmm (hmmpress'd)")
    s.add_argument("--threads", type=int, default=8)
    s.add_argument("--species", help="run one species only (for SLURM arrays)")
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-rescue", action="store_true", dest="no_rescue",
                   help="skip the below-threshold rescue pass for rescue=yes families")
    s.add_argument("--rescue-evalue", type=float, default=1e-3, dest="rescue_evalue")
    s.add_argument("--full-pfam", dest="full_pfam",
                   help="full Pfam-A.hmm for an architecture pass on candidate proteins")
    s.set_defaults(func=cmd_scan)

    n = sub.add_parser("annotate", help="apply domain rules, build count matrices")
    common(n)
    n.set_defaults(func=cmd_annotate)

    r = sub.add_parser("report", help="summaries and figures")
    common(r)
    r.add_argument("--metadata", help="TSV, first column or 'species' = proteome name")
    r.add_argument("--group-by", help="metadata column to compare groups by")
    r.set_defaults(func=cmd_report)

    e = sub.add_parser("extract", help="write per-family protein FASTA")
    common(e, True)
    e.add_argument("--families", help="comma-separated subset, default all")
    e.set_defaults(func=cmd_extract)

    d = sub.add_parser("domains", help="per-protein domain architecture and completeness")
    common(d)
    d.set_defaults(func=cmd_domains)

    q = sub.add_parser("qc", help="annotation quality from low-copy control families")
    common(q)
    q.set_defaults(func=cmd_qc)

    y = sub.add_parser("cyp", help="assign a family to clans using labelled references")
    common(y)
    y.add_argument("--refs", required=True,
                   help="FASTA with clan=XXX in each header (e.g. FlyBase CYP refs)")
    y.add_argument("--family", default="CYP", help="which family to classify")
    y.add_argument("--threads", type=int, default=8)
    y.set_defaults(func=cmd_cyp)

    t = sub.add_parser("trees", help="MAFFT + FastTree per family")
    common(t)
    t.add_argument("--families", help="comma-separated subset, default all extracted")
    t.add_argument("--max-seqs", type=int, default=6000, dest="max_seqs")
    t.add_argument("--threads", type=int, default=8)
    t.add_argument("--force", action="store_true")
    t.set_defaults(func=cmd_trees)

    w = sub.add_parser("tree", help="radial + rectangular species tree with defensome rings")
    common(w)
    w.add_argument("--tree", required=True, help="Newick species tree (e.g. OrthoFinder SpeciesTree_rooted.txt)")
    w.add_argument("--metadata"); w.add_argument("--group-by", dest="group_by")
    w.add_argument("--families", help="comma-separated ring families, default all CORE")
    w.add_argument("--phylogram", action="store_true", help="use branch lengths instead of a cladogram")
    w.add_argument("--raw", action="store_true", help="plot raw counts instead of per-10k")
    w.set_defaults(func=cmd_tree)

    gt = sub.add_parser("genetree", help="radial gene tree for one family")
    common(gt)
    gt.add_argument("--family", default="CYP")
    gt.add_argument("--tree", help="default results/trees/<FAMILY>.tre")
    gt.add_argument("--by", default="trait", choices=["trait","clan"])
    gt.add_argument("--samplesheet", help="colour by method from a sample sheet")
    gt.add_argument("--metadata"); gt.add_argument("--group-by", dest="group_by")
    gt.add_argument("--max-labels", type=int, default=400, dest="max_labels",
                    help="species labels only when the tree has at most this many tips")
    gt.set_defaults(func=cmd_genetree)

    cr = sub.add_parser("cyp-refs", help="build a clan-labelled CYP reference FASTA")
    cr.add_argument("--source", choices=["icpd", "flybase"], default="icpd")
    cr.add_argument("--fasta", help="protein FASTA (ICPD, or the FlyBase clan ref)")
    cr.add_argument("--table", help="ICPD sequence information table (xlsx/tsv/csv)")
    cr.add_argument("--orders", default="", help="e.g. Lepidoptera")
    cr.add_argument("--species", default="")
    cr.add_argument("--min-len", type=int, default=400, dest="min_len")
    cr.add_argument("--max-per-clan", type=int, default=300, dest="max_per_clan")
    cr.add_argument("--id-col", dest="id_col",
                    help="table column holding the FASTA header IDs; auto-detected by matching")
    cr.add_argument("--include-fragments", action="store_true", dest="include_fragments",
                    help="keep entries marked fragmented (default: exclude them)")
    cr.add_argument("--out", required=True)
    cr.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cr.add_argument("--proteomes", help=argparse.SUPPRESS)
    cr.set_defaults(func=cmd_cyp_refs, out_is_dir=False)

    ct = sub.add_parser("chem-triage", help="triage P450RDB before committing to it")
    ct.add_argument("--dir", required=True, help="directory holding the P450RDB CSVs")
    ct.add_argument("--out", required=True)
    ct.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    ct.add_argument("--proteomes", help=argparse.SUPPRESS)
    ct.set_defaults(func=cmd_chem_triage)

    ie = sub.add_parser("icpd-evidence",
                        help="sort ICPD reference records into evidence tiers")
    ie.add_argument("--tables", nargs="+", required=True,
                    help="References1.xlsx References2.xlsx")
    ie.add_argument("--out", required=True)
    ie.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    ie.add_argument("--proteomes", help=argparse.SUPPRESS)
    ie.set_defaults(func=cmd_icpd_evidence)

    cb = sub.add_parser("cyp-benchmark",
                        help="measure clan-assignment accuracy on held-out ICPD sequences")
    cb.add_argument("--fasta", required=True, help="ICPD protein FASTA")
    cb.add_argument("--table", required=True, help="ICPD information table")
    cb.add_argument("--clan-hmms", required=True, dest="clan_hmms",
                    help="directory of clan HMMs, e.g. results/cyp/clan_hmms")
    cb.add_argument("--refs", help="the reference FASTA, so its sequences are excluded")
    cb.add_argument("--orders", default="", help="restrict the test set, e.g. Lepidoptera")
    cb.add_argument("--min-len", type=int, default=400, dest="min_len")
    cb.add_argument("--threads", type=int, default=8)
    cb.add_argument("--out", required=True)
    cb.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cb.add_argument("--proteomes", help=argparse.SUPPRESS)
    cb.set_defaults(func=cmd_cyp_benchmark)

    su = sub.add_parser("setup", help="extract the defensome HMMs from Pfam-A; record the release")
    su.add_argument("--pfam", help="existing Pfam-A.hmm or Pfam-A.hmm.gz")
    su.add_argument("--download", action="store_true", help="fetch Pfam-A from EBI")
    su.add_argument("--url", help="override the Pfam-A download URL")
    su.add_argument("--keep-full", action="store_true", dest="keep_full",
                    help="with --download, also keep an uncompressed Pfam-A for --full-pfam")
    su.add_argument("--db-dir", default="db", dest="db_dir")
    su.add_argument("--force", action="store_true")
    su.add_argument("--map", default=DEFAULT_MAP)
    su.add_argument("--out", default=".", help=argparse.SUPPRESS)
    su.add_argument("--proteomes", help=argparse.SUPPRESS)
    su.set_defaults(func=cmd_setup, out_is_dir=False)

    dr = sub.add_parser("doctor", help="check tools and libraries; report which steps can run")
    dr.add_argument("--require", default="",
                    help="comma list (core,scan,cyp,trees); exit non-zero if any cannot run")
    dr.add_argument("--out", default=".", help=argparse.SUPPRESS)
    dr.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    dr.add_argument("--proteomes", help=argparse.SUPPRESS)
    dr.set_defaults(func=cmd_doctor, out_is_dir=False)

    insp = sub.add_parser("inspect",
                          help="examine proteome headers, write a draft sample sheet")
    insp.add_argument("--proteomes", required=True)
    insp.add_argument("--out", default=".")
    insp.add_argument("--samplesheet", help="output path, default <out>/samplesheet.tsv")
    insp.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    insp.set_defaults(func=cmd_inspect)

    col = sub.add_parser("collapse",
                         help="longest protein per gene, per the sample sheet")
    col.add_argument("--samplesheet", required=True)
    col.add_argument("--out", required=True)
    col.add_argument("--drop-utrorf", action="store_true", dest="drop_utrorf",
                     help="discard EvidentialGene utrorf entries")
    col.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    col.add_argument("--proteomes", help=argparse.SUPPRESS)
    col.set_defaults(func=cmd_collapse)

    cmp_ = sub.add_parser("compare",
                          help="compare defensome recovery across annotation methods")
    cmp_.add_argument("--samplesheet", required=True)
    cmp_.add_argument("--out", required=True)
    cmp_.add_argument("--reference", default="DToL",
                      help="method to express recovery against (default DToL)")
    cmp_.add_argument("--per10k", action="store_true",
                      help="normalise per 10k genes (default: raw, correct for paired within-species comparison)")
    cmp_.add_argument("--complete-only", action="store_true", dest="complete_only",
                      help="count only COMPLETE domain architectures; strips Trinity fragments")
    cmp_.add_argument("--exclude", default="",
                      help="comma-separated sample_ids or species to leave out")
    cmp_.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    cmp_.add_argument("--proteomes", help=argparse.SUPPRESS)
    cmp_.set_defaults(func=cmd_compare)

    so = sub.add_parser("short-orfs",
                        help="recover proteins too short for the ORF caller (metallothioneins)")
    so.add_argument("--samplesheet", required=True)
    so.add_argument("--out", required=True)
    so.add_argument("--nucleotide-col", default="nucleotide", dest="nucleotide_col",
                    help="sample-sheet column holding transcript FASTA paths")
    so.add_argument("--min-aa", type=int, default=25, dest="min_aa")
    so.add_argument("--max-aa", type=int, default=99, dest="max_aa",
                    help="upper bound; default 99, the band TransDecoder's default -m 100 discards")
    so.add_argument("--max-mb", type=float, default=2000, dest="max_mb")
    so.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    so.add_argument("--proteomes", help=argparse.SUPPRESS)
    so.set_defaults(func=cmd_short_orfs)

    dh = sub.add_parser("dataset-help", help="how to obtain and use ICPD, P450RDB and the human CYP dataset")
    dh.add_argument("--out", default=".", help=argparse.SUPPRESS)
    dh.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    dh.add_argument("--proteomes", help=argparse.SUPPRESS)
    dh.set_defaults(func=cmd_dataset_help)

    asx = sub.add_parser("assets", help="write the embedded map/dashboard files out for editing")
    asx.add_argument("--write", default=".", help="directory to write into")
    asx.add_argument("--force", action="store_true")
    asx.add_argument("--out", default=".", help=argparse.SUPPRESS)
    asx.add_argument("--map", default=DEFAULT_MAP, help=argparse.SUPPRESS)
    asx.add_argument("--proteomes", help=argparse.SUPPRESS)
    asx.set_defaults(func=cmd_assets)

    db = sub.add_parser("dashboard", help="single self-contained HTML dashboard")
    common(db)
    db.add_argument("--metadata"); db.add_argument("--group-by", dest="group_by")
    db.add_argument("--species-tree", dest="species_tree",
                    help="Newick species tree to embed (OrthoFinder SpeciesTree_rooted.txt)")
    db.add_argument("--dashboard-out", dest="dashboard_out",
                    help="output path, default <out>/dashboard.html")
    db.add_argument("--samplesheet",
                    help="sample sheet; adds species/method/pipeline, and makes method the default trait")
    db.add_argument("--light", action="store_true",
                    help="omit per-protein tables to keep the file small")
    db.set_defaults(func=cmd_dashboard)

    x = sub.add_parser("all", help="scan + annotate + report")
    common(x, True)
    x.add_argument("--pfam", required=True)
    x.add_argument("--threads", type=int, default=8)
    x.add_argument("--metadata")
    x.add_argument("--group-by")
    x.add_argument("--force", action="store_true")
    x.add_argument("--species", default=None)
    x.add_argument("--no-rescue", action="store_true", dest="no_rescue")
    x.add_argument("--rescue-evalue", type=float, default=1e-3, dest="rescue_evalue")
    x.add_argument("--full-pfam", dest="full_pfam")
    x.set_defaults(func=lambda a: (cmd_scan(a), cmd_annotate(a), cmd_report(a)))

    a = p.parse_args()
    # Most subcommands treat --out as a results directory; cyp-refs writes a
    # single FASTA there. Creating a directory with the file's name silently
    # breaks the write, so only mkdir when --out really is a directory.
    if getattr(a, "out_is_dir", True) and getattr(a, "out", None):
        os.makedirs(a.out, exist_ok=True)
    for opt in ("metadata", "group_by", "proteomes", "species", "families",
                "threads", "force", "refs", "family", "max_seqs", "tree",
                "phylogram", "raw", "by", "max_labels", "species_tree",
                "dashboard_out", "light", "write", "force", "source",
                "table", "orders", "min_len", "max_per_clan", "dir",
                "id_col", "include_fragments", "tables", "clan_hmms",
                "samplesheet", "drop_utrorf", "reference", "per10k",
                "complete_only", "exclude", "require", "no_rescue",
                "nucleotide_col", "min_aa", "max_aa", "max_mb",
                "rescue_evalue", "full_pfam", "download", "url", "keep_full",
                "db_dir"):
        if not hasattr(a, opt):
            setattr(a, opt, None)
    a.func(a)


if __name__ == "__main__":
    main()
