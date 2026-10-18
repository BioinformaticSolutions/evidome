"""What each database contributes, and how much the whole analysis explains.

The combination of databases is chosen by letting their gene sets compete, with
every order of entry considered (which database enters first matters), and a
goodness-of-fit score -- how much of the variation the finished analysis
explains -- shown as a curve.

Per comparison:

response    each tested human gene's own statistic -- the GSEA ranking metric
            (T-statistic, or sign x -log10 p) where the comparison has a
            direction, else -log10 p -- standardised
features    one column per gene set (Reactome, GO BP/MF/CC, WikiPathways,
            Hallmark; membership 0/1, 10-500 tested genes) and per regulator
            (CollecTRI, signalling via TF; the signed target vector, +1/-1,
            or membership where the response has no direction; 5+ tested
            targets, signalling at most 1000)
model       lasso regression of the response on the features of a set of
            databases; every set of every included database competes at once,
            so a set that repeats another explains nothing new and gets a zero
            coefficient
score v(S)  held-out R^2: 5-fold cross-validation over genes, at the penalty
            that predicts held-out genes best. Held out, a database of
            redundant or uninformative sets cannot raise the score -- the bias
            of in-sample R^2, which only grows as columns are added.
screening   inside each training fold, each database's columns are cut to the
            --screen (300) most correlated with the response (sure independence
            screening, Fan & Lv 2008), so every database enters with the same
            number of candidates and the fits stay fast; Hallmark's 50 all pass

The eight databases give 2^8 = 256 combinations; v(S) is computed for every
one. From them:

shapley     each database's average gain in v when it is added, over EVERY
            order of entry (exact Shapley value, Shapley 1953; as a variance
            decomposition, Lindeman, Merenda & Gold 1980) -- no order is
            privileged, so no sampling of orders is needed. The eight values
            add up to v(all).
goodness of fit
            v(all): the share of gene-to-gene variation in the comparison's
            statistic that the whole analysis (all databases together)
            predicts in genes it was not fitted on; with in-sample R^2, the
            number of sets and regulators the full model keeps, and a
            permutation null (the response shuffled across genes, --permutations
            20 times): its mean, 95th percentile and the empirical p
curve       rarefaction-style: for k = 1..8 databases, v over every
            combination of k (mean, lowest, highest) -- how much more is
            explained as databases are added, and where it levels off
stability   the full model refitted on --bootstrap (50) bootstrap samples of
            genes: how often each set or regulator is kept -- the redundancy-
            free list of what the analysis rests on

Caveat, recorded: genes sharing sets are not independent, so held-out R^2 is
somewhat optimistic in absolute terms; every database is scored the same way,
so comparisons between them are fair. This measures how much expression
change each database explains in this study, not whether its biology is true.

Requires numpy, scipy and scikit-learn (scikit-learn is a dependency and brings
the other two); should they be missing, the command says so (EVD-0028).
"""

from __future__ import annotations

import math
import os
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError

PATHWAY_DATABASES = ("reactome", "go_bp", "go_mf", "go_cc", "wikipathways", "hallmark")
REGULATOR_DATABASES = ("collectri", "signalling_via_tf")
DATABASES = PATHWAY_DATABASES + REGULATOR_DATABASES
SHAPLEY_COLUMNS = ("comparison", "database", "shapley_r2", "share_of_explained", "alone_r2", "loss_if_dropped",
                   "features")
FIT_COLUMNS = ("comparison", "response", "genes", "heldout_r2_all", "heldout_r2_pathways", "heldout_r2_regulators",
               "insample_r2_all", "features_kept", "null_mean_r2", "null_95th_r2", "null_p")
CURVE_COLUMNS = ("comparison", "databases", "combinations", "mean_r2", "lowest_r2", "highest_r2", "best_combination")
STABLE_COLUMNS = ("comparison", "feature", "database", "name", "selection_frequency", "coefficient")


@dataclass
class ValueSettings:
    folds: int = 5
    alphas: int = 12
    screen: int = 300
    bootstrap: int = 50
    permutations: int = 20
    seed: int = 1
    min_size: int = 10
    max_size: int = 500
    min_targets: int = 5
    max_targets: int = 1000


def run_limits(crosswalk_dir: Path) -> Dict[str, int]:
    """The gene-set size and regulator target limits the run itself used, from enrichment.json and
    regulators.json, so database-value scores what the run used rather than the defaults. Missing
    manifests leave the defaults."""
    import json

    out: Dict[str, int] = {}
    for rel, keys in (("enrichment/enrichment.json", (("min_size", "min_size"), ("max_size", "max_size"))),
                      ("regulators/regulators.json", (("min_targets", "min_targets"),
                                                      ("max_targets", "max_targets_signalling_via_tf")))):
        try:
            settings = json.loads((Path(crosswalk_dir) / rel).read_text(encoding="utf-8")).get("settings", {})
        except (OSError, ValueError, AttributeError):
            continue
        for mine, theirs in keys:
            if isinstance(settings.get(theirs), int):
                out[mine] = settings[theirs]
    return out


#: Maths libraries start one thread per core in EVERY process unless told otherwise; with one worker process per
#: core that is cores x cores threads fighting for the machine: three comparisons have run
#: 2 h 20 min on 40-core nodes without finishing one (each needs about 15 CPU-minutes). Set before numpy loads,
#: in the command and in the SLURM script, and workers are started fresh (spawn), never forked from a process
#: whose thread pools are already running.
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                 "VECLIB_MAXIMUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}


def single_thread_maths() -> None:
    for key, value in SINGLE_THREAD.items():
        os.environ[key] = value


def available() -> bool:
    try:
        import numpy  # noqa: F401
        import scipy  # noqa: F401
        import sklearn  # noqa: F401
    except ImportError:
        return False
    return True


def require() -> None:
    if not available():
        raise ProjectError(28, "numpy, scipy and scikit-learn are not all installed for this Python")


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def response(table) -> Tuple[List[str], List[float], str]:
    """(genes, values, description) for one comparison (enrich.ComparisonTable)."""
    genes = list(table.genes)
    if table.metric:
        return genes, [float(table.metric.get(g, 0.0)) for g in genes], table.metric_name
    return genes, [-math.log10(max(table.p[g], 1e-300)) for g in genes], "-log10(p) (no direction)"


def feature_blocks(genes: Sequence[str], info, members, regulons, settings: ValueSettings, signed: bool):
    """{database: (csc matrix genes x features, [feature ids], [names])}; databases with no feature are left out."""
    import numpy as np
    from scipy import sparse

    index = {g: i for i, g in enumerate(genes)}
    cols: Dict[str, Tuple[List[int], List[int], List[float], List[str], List[str]]] = {}
    for sid in sorted(members):
        db = str(info[sid]["source"])
        rows = [index[g] for g in members[sid] if g in index]
        if not settings.min_size <= len(rows) <= settings.max_size:
            continue
        r, c, v, ids, names = cols.setdefault(db, ([], [], [], [], []))
        c.extend([len(ids)] * len(rows))
        r.extend(rows)
        v.extend([1.0] * len(rows))
        ids.append(sid)
        names.append(str(info[sid].get("name", sid)))
    for (network, regulator), entry in sorted((regulons or {}).items()):
        targets = {g: s for g, s in entry["targets"].items() if g in index}
        if len(targets) < settings.min_targets or (network == "signalling_via_tf" and len(targets) > settings.max_targets):
            continue
        values = {g: float(s) for g, s in targets.items() if s} if signed else {g: 1.0 for g in targets}
        if not values:
            continue
        r, c, v, ids, names = cols.setdefault(network, ([], [], [], [], []))
        c.extend([len(ids)] * len(values))
        r.extend(index[g] for g in values)
        v.extend(values.values())
        ids.append(f"{network}:{regulator}")
        names.append(regulator)
    blocks = {}
    for db in DATABASES:
        if db in cols:
            r, c, v, ids, names = cols[db]
            blocks[db] = (sparse.csc_matrix((np.array(v), (np.array(r), np.array(c))), shape=(len(genes), len(ids))),
                          ids, names)
    return blocks


def _standardise(y):
    import numpy as np

    y = np.asarray(y, dtype=float)
    sd = y.std()
    return (y - y.mean()) / (sd if sd > 0 else 1.0)


def _screen(X, y, k: int):
    """Indices of the k columns of X most correlated (absolute) with y."""
    import numpy as np

    n = X.shape[0]
    yc = y - y.mean()
    cov = np.asarray(X.T @ yc).ravel() / n
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(sq - mean * mean, 1e-12))
    score = np.abs(cov) / sd
    if X.shape[1] <= k:
        return np.arange(X.shape[1])
    return np.sort(np.argpartition(-score, k - 1)[:k])


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

class Folds:
    """Cross-validation folds with each database's columns screened on the training genes."""

    def __init__(self, blocks, y, settings: ValueSettings, seed: int):
        import numpy as np
        from sklearn.model_selection import KFold

        self.y = y
        self.parts = []
        for train, test in KFold(settings.folds, shuffle=True, random_state=seed).split(y):
            per_db = {}
            for db, (X, _, _) in blocks.items():
                keep = _screen(X[train], y[train], settings.screen)
                sub = X[:, keep]
                per_db[db] = (sub[train].tocsc(), sub[test].tocsc())
            self.parts.append((train, test, per_db))
        top = max((float(np.max(np.abs(X.T @ (y - y.mean())))) / len(y) for X, _, _ in blocks.values()), default=1.0)
        self.alphas = np.geomspace(max(top, 1e-6), max(top, 1e-6) * 1e-3, settings.alphas)


def _fit_path(Xtr, ytr, alphas):
    """Coefficient paths from the largest penalty down (warm starts); yields (alpha index, model)."""
    import warnings

    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import Lasso

    model = Lasso(alpha=float(alphas[0]), fit_intercept=True, warm_start=True, max_iter=2000, tol=1e-3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        for i, a in enumerate(alphas):
            model.set_params(alpha=float(a))
            model.fit(Xtr, ytr)
            yield i, model


def heldout(folds: Folds, dbs: Sequence[str]) -> Tuple[float, int]:
    """v(S): the best mean held-out R^2 over the penalty grid, and that penalty's index. v(empty) = 0."""
    import numpy as np
    from scipy import sparse

    present = [d for d in dbs if d in folds.parts[0][2]]
    if not present:
        return 0.0, 0
    scores = np.zeros((len(folds.parts), len(folds.alphas)))
    for f, (train, test, per_db) in enumerate(folds.parts):
        Xtr = sparse.hstack([per_db[d][0] for d in present]).tocsc()
        Xte = sparse.hstack([per_db[d][1] for d in present]).tocsc()
        ytr, yte = folds.y[train], folds.y[test]
        sst = float(((yte - ytr.mean()) ** 2).sum())
        if sst == 0:            # nothing to explain in this fold: no credit
            continue
        for i, model in _fit_path(Xtr, ytr, folds.alphas):
            scores[f, i] = 1.0 - float(((yte - model.predict(Xte)) ** 2).sum()) / sst
    mean = scores.mean(axis=0)
    best = int(np.argmax(mean))
    # Floored at 0: predicting the training mean -- the empty model, always available -- scores 0 here, and the
    # penalty grid does not always reach it inside a fold; without the floor a useless database could score below
    # v(empty) and Shapley shares would pick up fitting noise.
    return max(0.0, float(mean[best])), best


def shapley(values: Dict[Tuple[str, ...], float], players: Sequence[str]) -> Dict[str, float]:
    """Exact Shapley values from v over every subset (keys: sorted tuples of players)."""
    n = len(players)
    out = {}
    for p in players:
        others = [q for q in players if q != p]
        total = 0.0
        for k in range(len(others) + 1):
            weight = math.factorial(k) * math.factorial(n - k - 1) / math.factorial(n)
            for S in combinations(others, k):
                total += weight * (values[tuple(sorted(S + (p,)))] - values[tuple(sorted(S))])
        out[p] = total
    return out


# workers: the fold data are sent once per process, then every subset is scored against them
_STATE: Dict[str, object] = {}


def _init(state) -> None:
    single_thread_maths()
    _STATE.update(state)


def _subset_job(dbs: Tuple[str, ...]) -> Tuple[Tuple[str, ...], float, int]:
    v, best = heldout(_STATE["folds"], dbs)  # type: ignore[arg-type]
    return dbs, v, best


def _null_job(seed: int) -> float:
    import numpy as np

    rng = np.random.default_rng(seed)
    y = rng.permutation(_STATE["y"])  # type: ignore[arg-type]
    folds = Folds(_STATE["blocks"], y, _STATE["settings"], seed)  # type: ignore[arg-type]
    return heldout(folds, tuple(_STATE["blocks"]))[0]  # type: ignore[arg-type]


def _boot_job(seed: int) -> List[str]:
    import numpy as np
    from scipy import sparse

    blocks, y, settings, alpha = _STATE["blocks"], _STATE["y"], _STATE["settings"], _STATE["alpha"]
    rng = np.random.default_rng(seed)
    rows = rng.integers(0, len(y), len(y))
    kept_ids: List[str] = []
    mats, ids = [], []
    for db, (X, fids, _) in blocks.items():  # type: ignore[union-attr]
        keep = _screen(X[rows], y[rows], settings.screen)  # type: ignore[union-attr]
        mats.append(X[rows][:, keep])
        ids.extend(fids[i] for i in keep)
    model = list(_fit_path(sparse.hstack(mats).tocsc(), y[rows], np.array([alpha])))[-1][1]
    kept_ids.extend(ids[i] for i in np.flatnonzero(model.coef_))
    return kept_ids


class _InProcess:
    """``--threads 1``: run each job here, in order, as parallel.map_ordered does."""

    def __init__(self, state):
        _init(state)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def submit(self, fn, *args):
        from concurrent.futures import Future

        future: Future = Future()
        future.set_result(fn(*args))
        return future

    def map(self, fn, items):
        return [fn(item) for item in items]


def _pool(threads: int, state):
    """The one exception to parallel.map_ordered (documented there): the workers need an initializer that
    loads the shared design matrices once, and the spawn context, so each starts with its maths libraries held
    to one thread. ``threads == 1`` runs in this process, as map_ordered does."""
    from concurrent.futures import ProcessPoolExecutor

    import multiprocessing

    if max(1, threads) == 1:
        return _InProcess(state)
    return ProcessPoolExecutor(max_workers=threads, initializer=_init, initargs=(state,),
                               mp_context=multiprocessing.get_context("spawn"))


def _progress(done: int, total: int, started: float, label: str, log) -> None:
    import time

    elapsed = time.time() - started
    left = elapsed / done * (total - done) if done else 0.0
    log(f"{label}: {done} of {total} done, {elapsed / 60:.1f} min so far, about {left / 60:.1f} min to go")


def value_comparison(name: str, kind: str, crosswalk_dir: Path, info, members, regulons, out: Path,
                     settings: ValueSettings, threads: int, log=print) -> Dict[str, object]:
    """Everything for one comparison, written under ``out/<name>/``. Returns its summary."""
    import numpy as np
    from scipy import sparse

    from .enrich import read_comparison

    table = read_comparison(Path(crosswalk_dir) / "comparisons" / f"{name}.human_genes.csv", name, kind)
    genes, raw, what = response(table)
    y = _standardise(raw)
    blocks = feature_blocks(genes, info, members, regulons, settings, signed=bool(table.metric))
    players = [d for d in DATABASES if d in blocks]
    log(f"{name}: {len(genes)} genes, response {what}; features "
        + ", ".join(f"{d} {blocks[d][0].shape[1]}" for d in players))
    if not players or len(set(raw)) < 2:
        # A response with no variation has nothing to explain.
        # No database contributed a single feature (every set outside the size range, no regulator index, or
        # --sources left nothing): nothing can be fitted, and the run says so.
        log(f"{name}: " + ("no database has a feature for these genes" if not players
                           else "the response does not vary") + "; nothing to fit")
        fit = {"comparison": name, "response": what, "genes": len(genes), "heldout_r2_all": 0.0,
               "heldout_r2_pathways": 0.0, "heldout_r2_regulators": 0.0, "insample_r2_all": 0.0,
               "features_kept": 0, "null_mean_r2": "", "null_95th_r2": "", "null_p": ""}
        folder = Path(out) / name
        folder.mkdir(parents=True, exist_ok=True)
        _write(folder / "shapley.csv", SHAPLEY_COLUMNS, [])
        _write(folder / "curve.csv", CURVE_COLUMNS, [])
        _write(folder / "fit.csv", FIT_COLUMNS, [fit])
        _write(folder / "stable_features.csv", STABLE_COLUMNS, [])
        _write(folder / "combinations.csv", ("comparison", "combination", "databases", "heldout_r2"), [])
        return {"fit": fit, "shapley": [], "curve": []}
    folds = Folds(blocks, y, settings, settings.seed)
    subsets = [tuple(sorted(S)) for k in range(1, len(players) + 1) for S in combinations(players, k)]
    state = {"folds": folds, "blocks": blocks, "y": y, "settings": settings}
    values: Dict[Tuple[str, ...], float] = {(): 0.0}
    best_alpha: Dict[Tuple[str, ...], int] = {}
    import time
    from concurrent.futures import as_completed

    started = time.time()
    step = max(1, len(subsets) // 10)
    log(f"{name}: scoring {len(subsets)} combinations of databases with {threads} worker(s)")
    with _pool(threads, state) as pool:
        futures = [pool.submit(_subset_job, S) for S in subsets]
        for done, fut in enumerate(as_completed(futures), 1):
            dbs, v, best = fut.result()
            values[dbs] = v
            best_alpha[dbs] = best
            if done % step == 0 or done == len(subsets):
                _progress(done, len(subsets), started, f"{name} combinations", log)
    full = tuple(sorted(players))
    log(f"{name}: {len(subsets)} combinations scored; all databases held-out R^2 {values[full]:.4f}")

    alpha = float(folds.alphas[best_alpha[full]])
    # The full model is screened exactly as the held-out fits and the bootstraps are: its penalty was chosen on
    # screened matrices, so fitting it unscreened would give in-sample R^2, features_kept and coefficients from a
    # model that was never tuned. It is refitted at its chosen penalty.
    mats, ids = [], []
    for d in players:
        X, fids, _ = blocks[d]
        keep = _screen(X, y, settings.screen)
        mats.append(X[:, keep])
        ids.extend(fids[i] for i in keep)
    design = sparse.hstack(mats).tocsc()
    model = list(_fit_path(design, y, np.array([alpha])))[-1][1]
    insample = float(model.score(design, y))
    kept = int(np.count_nonzero(model.coef_))
    names = {i: n for d in players for i, n in zip(blocks[d][1], blocks[d][2])}
    dbs_of = {i: d for d in players for i in blocks[d][1]}
    coef = dict(zip(ids, model.coef_))

    state.update({"alpha": alpha})
    log(f"{name}: {settings.permutations} shuffled-data fits and {settings.bootstrap} bootstrap fits")
    with _pool(threads, state) as pool:
        null = list(pool.map(_null_job, [settings.seed + 1000 + i for i in range(settings.permutations)]))
        boots = list(pool.map(_boot_job, [settings.seed + 5000 + i for i in range(settings.bootstrap)]))
    counts: Dict[str, int] = {}
    for kept_ids in boots:
        for i in kept_ids:
            counts[i] = counts.get(i, 0) + 1

    phi = shapley(values, players)
    total = values[full]
    shapley_rows = [{
        "comparison": name, "database": d, "shapley_r2": round(phi[d], 6),
        "share_of_explained": round(phi[d] / total, 4) if total > 0 else "",
        "alone_r2": round(values[(d,)], 6),
        "loss_if_dropped": round(total - values[tuple(sorted(set(full) - {d}))], 6),
        "features": blocks[d][0].shape[1]} for d in players]
    curve_rows = []
    for k in range(1, len(players) + 1):
        # In subset order, with ties broken by that order: `values` fills in completion order, so without this
        # equal R^2 would give a different best_combination from run to run.
        vs = {S: values[S] for S in subsets if len(S) == k}
        best = max(vs, key=lambda S: (vs[S], -subsets.index(S)))
        curve_rows.append({"comparison": name, "databases": k, "combinations": len(vs),
                           "mean_r2": round(sum(vs.values()) / len(vs), 6), "lowest_r2": round(min(vs.values()), 6),
                           "highest_r2": round(max(vs.values()), 6), "best_combination": "+".join(best)})
    fam = lambda group: values.get(tuple(sorted(d for d in players if d in group)), 0.0)   # noqa: E731
    null_sorted = sorted(null)
    fit = {"comparison": name, "response": what, "genes": len(genes), "heldout_r2_all": round(total, 6),
           "heldout_r2_pathways": round(fam(PATHWAY_DATABASES), 6),
           "heldout_r2_regulators": round(fam(REGULATOR_DATABASES), 6),
           "insample_r2_all": round(insample, 6), "features_kept": kept,
           "null_mean_r2": round(sum(null) / len(null), 6) if null else "",
           "null_95th_r2": round(null_sorted[int(0.95 * (len(null) - 1))], 6) if null else "",
           "null_p": round((sum(1 for v in null if v >= total) + 1) / (len(null) + 1), 4) if null else ""}
    stable = sorted(({"comparison": name, "feature": i, "database": dbs_of[i], "name": names[i],
                      "selection_frequency": round(c / max(1, settings.bootstrap), 3),
                      "coefficient": round(float(coef.get(i, 0.0)), 6)} for i, c in counts.items()),
                    key=lambda r: (-r["selection_frequency"], -abs(r["coefficient"]), r["feature"]))

    folder = Path(out) / name
    folder.mkdir(parents=True, exist_ok=True)
    _write(folder / "shapley.csv", SHAPLEY_COLUMNS, shapley_rows)
    _write(folder / "curve.csv", CURVE_COLUMNS, curve_rows)
    _write(folder / "fit.csv", FIT_COLUMNS, [fit])
    _write(folder / "stable_features.csv", STABLE_COLUMNS, stable)
    _write(folder / "combinations.csv", ("comparison", "combination", "databases", "heldout_r2"),
           [{"comparison": name, "combination": "+".join(S), "databases": len(S), "heldout_r2": round(v, 6)}
            for S, v in sorted(values.items(), key=lambda kv: (len(kv[0]), kv[0])) if S])
    log(f"{name}: done in {(time.time() - started) / 60:.1f} min; held-out R^2 {total:.4f}; wrote {folder}")
    return {"fit": fit, "shapley": shapley_rows, "curve": curve_rows}


def _write(path: Path, columns, rows) -> Path:
    import csv

    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)
    return path


def _read(path: Path) -> List[Dict[str, str]]:
    import csv

    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


# ---------------------------------------------------------------------------
# a run: every comparison, then the pooled table and the figures
# ---------------------------------------------------------------------------

def _inputs(reference_dir: Path, release: int, sources, networks, log):
    from .genesets import load_gene_sets
    from .regulators import load_regulons

    info, members = load_gene_sets(reference_dir, release, sources)
    try:
        regulons = load_regulons(reference_dir, release, networks)
    except ProjectError:
        regulons = {}
        log("  no regulator index: regulators are left out")
    return info, members, regulons


def finalize(crosswalk_dir: Path, out: Path, names: Sequence[str], log=print, figure_settings=None) -> List[str]:
    """Pool the per-comparison tables into database_value.csv and draw the curves."""
    out = Path(out)
    shap, fits, curves = [], [], []
    for n in names:
        folder = out / n
        if not (folder / "fit.csv").exists():
            raise ProjectError(22, f"database value for {n!r} did not finish (no {folder / 'fit.csv'})")
        shap += _read(folder / "shapley.csv")
        fits += _read(folder / "fit.csv")
        curves += _read(folder / "curve.csv")
    _write(out / "database_value.csv", SHAPLEY_COLUMNS, shap)
    _write(out / "goodness_of_fit.csv", FIT_COLUMNS, fits)
    _write(out / "curves.csv", CURVE_COLUMNS, curves)
    lines = value_lines(shap, fits)
    # beside the run's other figures by default; under --out when one is given, so --out's figures never
    # overwrite the default run's
    default = (Path(crosswalk_dir) / "database_value").resolve()
    fig_dir = Path(crosswalk_dir) / "figures" if out.resolve() == default else out / "figures"
    from . import figures

    if curves and not figures.available():
        lines.append("figures skipped: matplotlib is not installed (python3 -m pip install matplotlib)")
    drawn = draw_curves(curves, shap, fits, fig_dir, figure_settings)
    if drawn:
        lines.append(f"drew {drawn} database-value figure(s) in {fig_dir}")
    lines.append(f"wrote {out / 'database_value.csv'}, {out / 'goodness_of_fit.csv'}, {out / 'curves.csv'}")
    return lines


def value_lines(shap: List[Dict[str, str]], fits: List[Dict[str, str]]) -> List[str]:
    lines = ["goodness of fit: share of gene-to-gene variation explained in held-out genes (all databases)"]
    for f in fits:
        lines.append(f"  {f['comparison']:<34} R2 {float(f['heldout_r2_all']):.3f}  (pathways "
                     f"{float(f['heldout_r2_pathways']):.3f}, regulators {float(f['heldout_r2_regulators']):.3f}; "
                     f"null 95th {f['null_95th_r2']}, p {f['null_p']}; {f['features_kept']} features kept)")
    lines.append("database value: Shapley share of the explained variation, mean over comparisons")
    by: Dict[str, List[float]] = {}
    for r in shap:
        if r["share_of_explained"] != "":
            by.setdefault(r["database"], []).append(float(r["share_of_explained"]))
    for d in DATABASES:
        if d in by:
            lines.append(f"  {d:<19} {100 * sum(by[d]) / len(by[d]):6.1f}%")
    return lines


def draw_curves(curves, shap, fits, folder: Path, figure_settings=None) -> int:
    """Per comparison: the database curve (held-out R^2 against number of databases) beside the Shapley shares."""
    from . import figures

    if not figures.available() or not curves:
        return 0
    import matplotlib
    from matplotlib.figure import Figure

    settings = figure_settings or figures.FigureSettings()
    low, high = figures.parse_colors(settings.colors)
    folder.mkdir(parents=True, exist_ok=True)
    names = list(dict.fromkeys(r["comparison"] for r in curves))
    rc = {**figures.REPRODUCIBLE_RC, "font.family": "sans-serif", "font.sans-serif": figures.FONT, "font.size": 7, "axes.linewidth": 0.5,
          "xtick.major.width": 0.5, "ytick.major.width": 0.5, "pdf.fonttype": 42, "svg.fonttype": "none"}
    drawn = 0
    with matplotlib.rc_context(rc):
        for name in names:
            cv = [r for r in curves if r["comparison"] == name]
            sh = sorted((r for r in shap if r["comparison"] == name), key=lambda r: float(r["shapley_r2"]))
            fit = next(f for f in fits if f["comparison"] == name)
            fig = Figure(figsize=(7.0, 2.8))
            ax, bx = fig.add_axes([0.08, 0.17, 0.40, 0.68]), fig.add_axes([0.66, 0.17, 0.31, 0.68])
            k = [int(r["databases"]) for r in cv]
            ax.fill_between(k, [float(r["lowest_r2"]) for r in cv], [float(r["highest_r2"]) for r in cv],
                            color=high, alpha=0.15, linewidth=0)
            ax.plot(k, [float(r["mean_r2"]) for r in cv], color=high, linewidth=1.2, marker="o", markersize=3)
            ax.plot(k, [float(r["highest_r2"]) for r in cv], color=high, linewidth=0.6, linestyle="--")
            if fit.get("null_95th_r2") not in ("", None):
                ax.axhline(float(fit["null_95th_r2"]), color="#808080", linewidth=0.6, linestyle=":")
            ax.set_xlabel("Databases included")
            ax.set_ylabel("Held-out R$^2$")
            ax.set_xticks(k)
            ax.set_title("Variation explained as databases are added", fontsize=7.5, pad=4)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ys = list(range(len(sh)))
            vals = [float(r["shapley_r2"]) for r in sh]
            bx.barh(ys, vals, color=[high if v >= 0 else low for v in vals], height=0.65)
            bx.set_yticks(ys)
            bx.set_yticklabels([r["database"] for r in sh])
            bx.axvline(0, color="black", linewidth=0.5)
            bx.set_xlabel("Shapley share of held-out R$^2$")
            bx.set_title("Contribution of each database", fontsize=7.5, pad=4)
            for s in ("top", "right"):
                bx.spines[s].set_visible(False)
            fig.text(0.5, 0.97, f"{name} — goodness of fit: held-out R$^2$ {float(fit['heldout_r2_all']):.3f} "
                                f"(null p {fit['null_p']})", ha="center", va="top", fontsize=8)
            for fmt in settings.formats:
                figures.save_reproducibly(fig, folder / f"{name}_database_value.{fmt}", fmt, settings.dpi)
            drawn += 1
    return drawn


def database_value(crosswalk_dir: Path, reference_dir: Path, *, release: int = 116, out: Optional[Path] = None,
                   settings: Optional[ValueSettings] = None, sources=None, networks=None,
                   threads: Optional[int] = None, figure_settings=None, log=print) -> List[str]:
    """Every comparison here, one after another (each uses every CPU), then finalize."""
    from .enrich import _comparisons
    from .parallel import resolve_threads

    require()
    settings = settings or ValueSettings()
    out = Path(out) if out else Path(crosswalk_dir) / "database_value"
    comps = _comparisons(Path(crosswalk_dir))
    info, members, regulons = _inputs(reference_dir, release, sources, networks, log)
    for name, kind in comps:
        value_comparison(name, kind, Path(crosswalk_dir), info, members, regulons, out, settings,
                         resolve_threads(threads), log)
    return finalize(crosswalk_dir, out, [n for n, _ in comps], log, figure_settings)


PLAN = "plan.json"


def plan_slurm(crosswalk_dir: Path, reference_dir: Path, *, options, release: int = 116, out: Optional[Path] = None,
               settings: Optional[ValueSettings] = None, sources=None, networks=None, figure_settings=None,
               submit=None, log=print) -> Dict[str, str]:
    """One array task per comparison, then a finalize job after them."""
    import shlex
    import sys

    from .enrich import _comparisons
    from .identity import NAME
    from .jsonio import write_json
    from .slurm import script_text, submit as sbatch_submit
    from .version import __version__

    require()
    submit = submit or sbatch_submit
    settings = settings or ValueSettings()
    crosswalk_dir, reference_dir = Path(crosswalk_dir).resolve(), Path(reference_dir).resolve()
    comps = _comparisons(crosswalk_dir)
    out = (Path(out) if out else crosswalk_dir / "database_value").resolve()
    work = out / "_slurm"
    (work / "logs").mkdir(parents=True, exist_ok=True)
    plan = write_json(work / PLAN, {
        "planned_by": __version__, "crosswalk_dir": str(crosswalk_dir), "reference_dir": str(reference_dir),
        "out": str(out), "release": release, "settings": asdict(settings), "sources": sources,
        "networks": networks, "figures": figure_settings.as_dict() if figure_settings else None,
        "comparisons": [{"name": n, "kind": k} for n, k in comps]})
    command = " ".join(f"{k}={v}" for k, v in SINGLE_THREAD.items()) + f" {shlex.quote(sys.executable)} -m {NAME}.cli"
    task = script_text(f"{NAME}-dbvalue", options, work / "logs" / "task.%A_%a.log",
                       f"{command} database-value-task --plan {shlex.quote(str(plan))} --index \"$SLURM_ARRAY_TASK_ID\"",
                       array=f"0-{len(comps) - 1}")
    final = script_text(f"{NAME}-dbvalue-finalize", options.small(), work / "logs" / "finalize.%j.log",
                        f"{command} database-value-finalize --plan {shlex.quote(str(plan))}")
    (work / "task.sbatch").write_text(task)
    (work / "finalize.sbatch").write_text(final)
    array_id = submit(work / "task.sbatch")
    final_id = submit(work / "finalize.sbatch", dependency=f"afterok:{array_id}")
    log(f"submitted {len(comps)} comparison task(s) as array job {array_id} ({options.cpus_per_task} CPUs each), "
        f"finalize job {final_id} after it")
    log(f"output and logs: {out}  (progress: tail -f {shlex.quote(str(work / 'logs'))}/task.{array_id}_*.log)")
    return {"array": array_id, "finalize": final_id, "plan": str(plan), "logs": str(work / "logs")}


def _read_plan(path: Path):
    from .jsonio import read_json
    from .version import __version__

    plan = read_json(Path(path))
    if plan["planned_by"] != __version__:
        raise ProjectError(22, f"this node runs {__version__} but the plan was written by {plan['planned_by']}: "
                               "install the same version on every node")
    return plan, ValueSettings(**plan["settings"])


def run_task(plan_path: Path, index: int, *, threads: Optional[int] = None, log=print) -> Dict[str, object]:
    from .parallel import resolve_threads

    require()
    plan, settings = _read_plan(plan_path)
    comps = plan["comparisons"]
    if not 0 <= index < len(comps):
        raise ProjectError(22, f"task index {index} outside 0..{len(comps) - 1}")
    info, members, regulons = _inputs(Path(plan["reference_dir"]), plan["release"], plan["sources"],
                                      plan["networks"], log)
    c = comps[index]
    return value_comparison(c["name"], c["kind"], Path(plan["crosswalk_dir"]), info, members, regulons,
                            Path(plan["out"]), settings, resolve_threads(threads), log)


def run_finalize(plan_path: Path, log=print) -> List[str]:
    from .figures import FigureSettings

    plan, _ = _read_plan(plan_path)
    fs = FigureSettings.from_dict(plan["figures"]) if plan.get("figures") else None
    return finalize(Path(plan["crosswalk_dir"]), Path(plan["out"]), [c["name"] for c in plan["comparisons"]], log, fs)

