"""Reproducible explanatory figures and illustrative model outputs.

The worked-example model below mirrors the C# implementation in
run_verification.ps1 (ExampleInp, Period, RunDesign);
check_against_script() asserts that it reproduces the values reported there.
"""
import csv
import math
from numbers import Real
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image

OUT = Path(".")
plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})

K = 6
iF, iD, iO, iL, iR, iT = range(6)
DESIGNS = {
    "A": dict(mcov=0.5, macc=0.8, mtime=0.5, mplur=0.3, gamma=0.2, rho=0.2),
    "B": dict(mcov=0.9, macc=0.8, mtime=0.9, mplur=0.8, gamma=0.9, rho=0.8),
}
BASE = dict(delta0=0.02, eta=0.6, h=2.0, theta=0.5, phi=0.5, chi=0.02, lockR=0.005, q=0.5, s0=0.4,
            kappa=0.25, sig=0.8, noise=0.2, c0=0.4, period=1.0, cap=1.5,
            routine=(0.0, 0.0), floors=(0.0, 0.0, 0.0, 0.0, 0.30, 0.0),
            monitoring_form="product", salience_form="hill", horizon_tolerance=0.01,
            required_routes=(True,), available_routes=(True,))
EPS_S, EPS_C, EPS_D, MQMIN = 0.01, 0.0, 0.05, 0.10
SIMULATION_EVIDENCE_STATUS = "NotAssessed_IllustrativeSimulation"
LATENT_CLASS_NAMES = dict(X="Unacceptable", U="Unverified", F="Fragile", G="CapabilityGain",
                         R="Robust", S="Resilient")


def clip(x):
    finite(x, "projection input")
    return max(0.0, min(1.0, x))

def finite(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number.")
    return value

def unit(value, name):
    finite(value, name)
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be within [0,1].")
    return value


def vector(values, size, name, *, normalized=False):
    if values is None or not hasattr(values, "__len__") or len(values) != size:
        raise ValueError(f"{name} must contain exactly {size} values.")
    for value in values:
        (unit if normalized else finite)(value, name)


def stable_norm(values):
    if values is None:
        raise ValueError("A finite vector is required.")
    values = tuple(values)
    for value in values:
        finite(value, "norm coordinate")
    result = math.hypot(*values)
    return finite(result, "norm result")


def salience(strength, scale, exponent):
    finite(strength, "evidence strength")
    finite(scale, "salience scale")
    finite(exponent, "Hill exponent")
    if strength < 0 or scale <= 0 or exponent < 1:
        raise ValueError("Require strength >= 0, scale > 0 and Hill exponent >= 1.")
    if strength == 0:
        return 0.0
    if 0.5 * scale <= strength <= 2 * scale:
        log_ratio = math.log1p((strength - scale) / scale)
    else:
        log_ratio = math.log(strength) - math.log(scale)
    z = exponent * log_ratio
    if z >= 0:
        return 1 / (1 + math.exp(-z))
    ez = math.exp(z)
    return ez / (1 + ez)


def review_rate(annual, interval):
    unit(annual, "annual rate")
    finite(interval, "review interval")
    if interval <= 0:
        raise ValueError("The review interval must be positive.")
    if annual == 0 or annual == 1 or interval == 1:
        return annual
    return -math.expm1(interval * math.log1p(-annual))


def _capacity_coordinate(stock, eta, retained, relevance, valid, accuracy, decay, consumption):
    learning = eta * retained * relevance
    invalidation = learning * (1 - valid) * (1 - accuracy)
    budget = eta + decay + consumption
    if budget > 1 + 4 * math.ulp(1.0):
        raise ArithmeticError("The capacity survival coefficient violates the admitted rate bound.")
    # S29 bounds the coefficient, not the stock; absorb only coefficient-summation roundoff.
    survival = (1 - min(1.0, budget)) + (eta - invalidation)
    return learning * valid * (1 - stock) + survival * stock


def validity(signal, error):
    if signal is None or error is None or len(signal) != len(error):
        raise ValueError("Signal and error must use the same canonical evidence space.")
    a, b = stable_norm(signal), stable_norm(error)
    scale = max(a, b)
    return (a / scale) / (a / scale + b / scale) if scale else 0.0


def risk_accepts(bounds, tolerances):
    if bounds is None or tolerances is None or len(bounds) != len(tolerances) or not bounds:
        raise ValueError("Every aggregate/group/joint bound requires a corresponding tolerance.")
    for bound, tolerance in zip(bounds, tolerances):
        unit(bound, "risk bound")
        unit(tolerance, "risk tolerance")
        if tolerance == 1:
            raise ValueError("Risk tolerances must lie in [0,1), not [0,1].")
    return all(bound <= tolerance for bound, tolerance in zip(bounds, tolerances))


def critical_learning(eta, qv, erosion):
    unit(eta, "eta")
    finite(qv, "Qv")
    finite(erosion, "erosion")
    if erosion < 0:
        raise ValueError("Erosion must be nonnegative.")
    if abs(qv) > 1:
        raise ValueError("Qv must be within [-1,1].")
    if eta == 0 or qv <= 0:
        return dict(value=None, status="NoPositiveGainAtAnyLearningLevel")
    denominator = eta * qv
    if denominator == 0:
        raise ValueError("Critical-learning denominator underflowed; higher precision is required.")
    threshold = erosion / denominator
    return dict(value=threshold if math.isfinite(threshold) else None,
                status=("UnattainableWithinUnitLearning" if denominator <= erosion
                        else "FiniteStrictThreshold"))

def invalid_lesson_sensitivity(rejected_invalid, retained_invalid):
    if any(type(value) is not int or value < 0 for value in (rejected_invalid, retained_invalid)):
        raise ValueError("Independently adjudicated lesson counts must be nonnegative integers.")
    total = rejected_invalid + retained_invalid
    return dict(value=rejected_invalid / total if total else None,
                status=("Assessed_InvalidLessonRejectionSensitivity" if total
                        else "Unassessed_NoIndependentlyInvalidLessons"))


def required_routes_hold(required, available):
    if required is None or available is None or len(required) != len(available):
        raise ValueError("Required reporting routes and explicit availability flags must agree.")
    if any(type(value) is not bool for value in (*required, *available)):
        raise ValueError("Reporting-route flags must be booleans, not inferred quality scores.")
    return all(not need or present for need, present in zip(required, available))


def canonical_evidence(values, references, channels, weights, source_ids):
    """Fuse one reconciled contribution per original source/channel with frozen weights."""
    if (any(value is None or not hasattr(value, "__len__")
            for value in (values, references, channels, weights, source_ids)) or not len(values)
            or not len(values) == len(references) == len(channels) == len(weights) == len(source_ids)):
        raise ValueError("Evidence, reference, channel, fixed source-weight and source/version vectors must agree.")
    if any(type(c) is not int or c < 0 for c in channels):
        raise ValueError("Canonical channel identifiers must be nonnegative integers.")
    if max(channels) >= len(values):
        raise ValueError("Canonical channel identifiers must be contiguous from zero.")
    if sorted(set(channels)) != list(range(max(channels) + 1)):
        raise ValueError("Canonical channel identifiers must be contiguous from zero.")
    for value in values:
        finite(value, "evidence observation")
    for value in references:
        finite(value, "evidence reference")
    for value in weights:
        finite(value, "source weight")
    if any(r <= 0 for r in references):
        raise ValueError("Every evidence reference must be finite and positive.")
    if any(w < 0 for w in weights):
        raise ValueError("Fusion weights must be finite and nonnegative.")
    unique = {}
    for value, reference, channel, weight, source in zip(values, references, channels, weights, source_ids):
        if not isinstance(source, str) or not source.strip():
            raise ValueError("Stable original-source/version IDs are required.")
        normalized = value / reference
        finite(normalized, "normalized evidence")
        key = channel, source
        observation = normalized, weight
        if key in unique and unique[key] != observation:
            raise ValueError("Conflicting versions or weights must be reconciled before source fusion.")
        unique[key] = observation
    result = []
    for channel in range(max(channels) + 1):
        sources = [observation for (c, _), observation in unique.items() if c == channel]
        scale = max(weight for _, weight in sources)
        if scale <= 0:
            raise ValueError("Every canonical channel requires positive total source weight.")
        total = math.fsum(weight / scale for _, weight in sources)
        fused = math.fsum((weight / scale) / total * x for x, weight in sources)
        result.append(finite(fused, "fused evidence"))
    return tuple(result)


def monitoring_quality(components, form):
    vector(components, 3, "monitoring components", normalized=True)
    product = math.prod(components)
    if form == "product":
        return product
    if form == "geometric":
        return product ** (1 / 3)
    if form == "minimum":
        return min(components)
    raise ValueError(f"Unknown monitoring form: {form}")


def protection_holds(performance, performance_floors, capacity, capacity_floors,
                     group_harm, group_limits, *, peak=(), peak_floors=(),
                     cumulative=(), cumulative_limits=()):
    pairs = ((performance, performance_floors, True), (capacity, capacity_floors, True),
             (peak, peak_floors, True), (group_harm, group_limits, False),
             (cumulative, cumulative_limits, False))
    for observed, limits, normalized in pairs:
        if (observed is None or limits is None or not hasattr(observed, "__len__")
                or not hasattr(limits, "__len__") or len(observed) != len(limits)):
            raise ValueError("Every protected outcome requires a matching declared limit.")
        vector(observed, len(observed), "protected values", normalized=normalized)
        vector(limits, len(limits), "declared limits", normalized=normalized)
        if any(x < 0 for x in (*observed, *limits)):
            raise ValueError("Protected values and limits must be finite and in their declared domain.")
    return (all(x >= limit for x, limit in zip(performance, performance_floors))
            and all(x >= limit for x, limit in zip(capacity, capacity_floors))
            and all(x >= limit for x, limit in zip(peak, peak_floors))
            and all(x <= limit for x, limit in zip(group_harm, group_limits))
            and all(x <= limit for x, limit in zip(cumulative, cumulative_limits)))


def validate_period_inputs(D, P, C, S, dscale):
    if not isinstance(D, dict) or not isinstance(P, dict):
        raise ValueError("Explicit governance and worked-example parameter dictionaries are required.")
    if set(BASE) - P.keys() or set(DESIGNS["B"]) - D.keys():
        raise ValueError("All documented parameters, protections and governance inputs must be supplied.")
    vector(C, K, "capacity", normalized=True)
    vector(S, 2, "performance", normalized=True)
    vector(P["routine"], 2, "routine signal")
    vector(P["floors"], K, "capacity floors", normalized=True)
    for name in DESIGNS["B"]:
        unit(D[name], name)
    for name in ("delta0", "eta", "theta", "phi", "chi", "lockR", "q", "c0", "horizon_tolerance"):
        unit(P[name], name)
    for name in ("h", "s0", "kappa", "sig", "noise", "period"):
        finite(P[name], name)
    finite(dscale, "disturbance multiplier")
    if P["h"] < 1 or P["s0"] <= 0 or P["kappa"] <= 0 or P["period"] <= 0 or P["sig"] < 0 or dscale < 0:
        raise ValueError("Require h >= 1, positive scales/period and nonnegative disturbance.")
    if (isinstance(P["cap"], bool) or not isinstance(P["cap"], Real)
            or math.isnan(P["cap"]) or P["cap"] <= 0):
        raise ValueError("Saturation cap must be positive; infinity selects no cap.")
    required_routes_hold(P["required_routes"], P["available_routes"])
    if len(P["required_routes"]) != 1 or len(P["available_routes"]) != 1:
        raise ValueError("This one-receptor specialization requires one explicit reporting-route policy.")
    delta0 = review_rate(P["delta0"], P["period"])
    lock = review_rate(P["lockR"], P["period"])
    if P["eta"] + delta0 + P["chi"] + lock > 1:
        raise ValueError("S29 must hold after conversion to the actual review interval.")
    return delta0, lock


def period(D, P, C, S, dscale):
    """The two-performance-coordinate specialization, not a calibrated field model."""
    delta0, lockR = validate_period_inputs(D, P, C, S, dscale)
    d = min(P["cap"], P["sig"] * P["period"]) * dscale
    finite(d, "review disturbance")
    sig = (d + P["routine"][0], P["routine"][1])
    xi = (0.0, P["noise"])
    e = (sig[0] + xi[0], sig[1] + xi[1])
    v = validity(sig, xi)
    ne = stable_norm(e)
    mod = lambda c: P["theta"] + (1 - P["theta"]) * c
    MQ = monitoring_quality((D["mcov"], D["mtime"], D["mplur"] * mod(C[iD])), P["monitoring_form"])
    routes_pass = required_routes_hold(P["required_routes"], P["available_routes"])
    if not routes_pass:
        MQ = 0.0
    if P["salience_form"] == "hill":
        sal = salience(ne, P["s0"], P["h"])
    elif P["salience_form"] == "linear_capped":
        sal = min(1.0, ne / P["s0"])
    else:
        raise ValueError(f"Unknown salience form: {P['salience_form']}")
    ell = MQ * D["gamma"] * sal
    rhoT = D["rho"] * mod(C[iL])
    ellR = rhoT * ell
    a = 1 + 0.5 * C[iT]
    us = [ell * a * (e[0] / ne if ne > 0 else 0), ell * a * (e[1] / ne if ne > 0 else 0)]
    bound = P["kappa"] * (C[iF] + C[iO] + C[iR])
    nus = stable_norm(us)
    if nus > bound:
        us = [(u / nus) * bound if bound > 0 else 0.0 for u in us]
    nu = min(1.0, stable_norm(us) / bound) if bound > 0 else 0.0
    loss = 0.05 * d
    raw_proxy = [S[0] - loss, S[1]]
    raw_performance = [S[0] - loss + 0.2 * us[0], S[1]]
    Smin, Sn = list(map(clip, raw_proxy)), list(map(clip, raw_performance))
    deltaT = delta0 * (1 - P["phi"] * rhoT)
    zeta = [0.0] * K
    zeta[iO] = P["chi"] * nu
    zeta[iR] = P["chi"] * nu + lockR
    acc = D["macc"]
    Cn = [_capacity_coordinate(C[k], P["eta"], ellR, P["q"], v, acc, deltaT, zeta[k])
          for k in range(K)]
    vector(Cn, K, "capacity output", normalized=True)
    SB = protection_holds(Sn, (0.45, 0.0), Cn, P["floors"], (0.3 * d,), (0.3,))
    SC, SCn, SCmin = 0.2 * S[0] + 0.8 * S[1], 0.2 * Sn[0] + 0.8 * Sn[1], 0.2 * Smin[0] + 0.8 * Smin[1]
    return dict(Cn=Cn, Sn=Sn, MQ=MQ, ell=ell, ellR=ellR, v=v, SB=SB, SC=SC, SCn=SCn, SCmin=SCmin,
                Smin=Smin, d=d, evidence_strength=ne, capacity_floors=tuple(P["floors"]),
                raw_performance=raw_performance, raw_impact_proxy=raw_proxy, raw_capacity=list(Cn),
                required_routes_pass=routes_pass, path_evidence_status="ProxyOnly_NotObservedPathSafety")


def latent_class_of(R, Rz):
    for name in ("MQ", "SC", "SCn", "SCmin"):
        unit(R[name], name)
    finite(R["d"], "disturbance")
    if R["d"] < 0 or type(R["SB"]) is not bool:
        raise ValueError("Declare nonnegative disturbance and an explicit safeguard flag.")
    for state in (R["Cn"], Rz["C"], Rz["Cn"]):
        vector(state, K, "classification capacity", normalized=True)
    idx = lambda C: sum(C) / K
    dC = idx(R["Cn"]) - idx(Rz["C"])
    dCd = idx(R["Cn"]) - idx(Rz["Cn"])
    if not R["SB"]:
        return "X"
    if R["MQ"] < MQMIN:
        return "U"
    if R["SCn"] < R["SC"] - EPS_S:
        return "F"
    if R["d"] > EPS_D and dCd > EPS_C and dC >= -EPS_C:
        return "G"
    return "R" if R["SCmin"] >= R["SC"] - EPS_S else "S"

def review_protection_holds(result, initial_index, tolerance):
    for name in ("MQ", "SC", "SCn", "SCmin"):
        unit(result[name], name)
    unit(initial_index, "horizon-start performance")
    unit(tolerance, "horizon tolerance")
    if type(result["SB"]) is not bool:
        raise ValueError("An explicit safeguard flag is required.")
    return (result["SB"] and result["MQ"] >= MQMIN
            and result["SCn"] >= result["SC"] - EPS_S
            and result["SCn"] >= initial_index - tolerance)


def run_design(name, P, n, dist, C0=None, S0=(0.55, 0.60)):
    if name not in DESIGNS or type(n) is not int or n <= 0 or not callable(dist):
        raise ValueError("Declare a known design, positive integer review count and disturbance schedule.")
    D = DESIGNS[name]
    C = list(C0) if C0 is not None else [P["c0"]] * K
    Ccf = list(C)
    S = list(S0)
    Scf = list(S0)
    validate_period_inputs(D, P, C, S, 0.0)
    baseline = 0.2 * S[0] + 0.8 * S[1]
    performance, performance_cf = [tuple(S)], [tuple(Scf)]
    idx, idxcf, classes = [sum(C) / K], [sum(C) / K], []
    components, components_cf = [tuple(C)], [tuple(Ccf)]
    protected = protection_holds(S, (0.45, 0.0), C, P["floors"], (), ())
    material = False
    for t in range(n):
        on = dist(t)
        if type(on) is not bool:
            raise ValueError("The disturbance schedule must return an explicit boolean for each review.")
        R = period(D, P, C, S, 1.0 if on else 0.0)
        Z = period(D, P, C, S, 0.0)
        Z["C"] = C
        classes.append(latent_class_of(R, Z))
        review_safe = review_protection_holds(R, baseline, P["horizon_tolerance"])
        protected = protected and review_safe
        material = material or R["d"] > EPS_D
        C, S = R["Cn"], R["Sn"]
        comparator = period(D, P, Ccf, Scf, 0.0)
        Ccf, Scf = comparator["Cn"], comparator["Sn"]
        idx.append(sum(C) / K)
        idxcf.append(sum(Ccf) / K)
        components.append(tuple(C))
        components_cf.append(tuple(Ccf))
        performance.append(tuple(S))
        performance_cf.append(tuple(Scf))
    return dict(idx=idx, idxcf=idxcf, classes="".join(classes), C=C,
                components=components, components_cf=components_cf,
                performance=performance, performance_cf=performance_cf,
                initial_performance_index=baseline, horizon_tolerance=P["horizon_tolerance"],
                protected_horizon=protected,
                retained_gain=(protected and material and idx[-1] - idxcf[-1] > EPS_C
                               and idx[-1] - idx[0] >= -EPS_C),
                evidence_status=SIMULATION_EVIDENCE_STATUS)


def export_model_records():
    rows = []
    for convention, cap in (("capped", 1.5), ("uncapped", math.inf)):
        P = dict(BASE, period=5.0, cap=cap)
        for name in DESIGNS:
            result = period(DESIGNS[name], P, [0.4] * K, [0.55, 0.60], 1.0)
            zero = period(DESIGNS[name], P, [0.4] * K, [0.55, 0.60], 0.0)
            zero["C"] = [0.4] * K
            rows.append(dict(convention=convention, design=name, disturbance=result["d"],
                             social_performance=result["Sn"][0], group_harm=0.3 * result["d"],
                             capacity_index=sum(result["Cn"]) / K,
                             reversibility=result["Cn"][iR], reversibility_floor=P["floors"][iR],
                             latent_class=LATENT_CLASS_NAMES[latent_class_of(result, zero)],
                             evidence_status=SIMULATION_EVIDENCE_STATUS))
    with (OUT / "Cap_Comparison.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    rows = []
    for name in DESIGNS:
        result = run_design(name, dict(BASE), 30, lambda t: t % 5 == 2)
        for t in range(31):
            row = dict(design=name, period=t, capacity_index=result["idx"][t],
                       no_disturbance_index=result["idxcf"][t],
                       latent_class=LATENT_CLASS_NAMES[result["classes"][t - 1]] if t else "NotClassified",
                       evidence_status=result["evidence_status"],
                       performance_index=0.2 * result["performance"][t][0] + 0.8 * result["performance"][t][1],
                       counterfactual_performance_index=0.2 * result["performance_cf"][t][0] + 0.8 * result["performance_cf"][t][1],
                       horizon_start_performance_index=result["initial_performance_index"],
                       horizon_tolerance=result["horizon_tolerance"],
                       protected_horizon=result["protected_horizon"], retained_gain=result["retained_gain"],
                       path_evidence_status="ProxyOnly_NotObservedPathSafety")
            for k, capacity in enumerate(("F", "D", "O", "L", "R", "T")):
                row[capacity] = result["components"][t][k]
                row[f"counterfactual_{capacity}"] = result["components_cf"][t][k]
                row[f"floor_{capacity}"] = BASE["floors"][k]
            rows.append(row)
    with (OUT / "Illustrative_Trajectories.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check_against_script():
    P = dict(BASE)
    a = run_design("A", P, 5, lambda t: True)
    b = run_design("B", P, 5, lambda t: True)
    assert f"{b['idx'][5]:.4f}" == "0.4811" and f"{b['idxcf'][5]:.4f}" == "0.3612", (b["idx"][5], b["idxcf"][5])
    assert f"{a['idx'][5]:.4f}" == "0.3630" and a["classes"] == "UUXXX" and b["classes"] == "GGGGG"
    ep = lambda t: t % 5 == 2
    hb = run_design("B", P, 30, ep)
    assert f"{hb['idx'][30]:.4f}" == "0.3373" and f"{hb['idxcf'][30]:.4f}" == "0.2154", (hb["idx"][30], hb["idxcf"][30])
    long = run_design("B", P, 4000, ep)
    hm = run_design("B", P, 30, ep, C0=long["C"])
    assert f"{hm['idx'][0]:.4f}" == "0.2761" and f"{hm['idx'][30]:.4f}" == "0.2761", (hm["idx"][0], hm["idx"][30])
    print("Python model reproduces the script values.")


# ---------------------------------------------------------------- drawing helpers
BLUE, BLUE_FILL, GREY_FILL, EDGE = "#0072B2", "#DFEBF5", "#F1F1F1", "#3A3A3A"


def box(ax, x, y, w, h, text, fill=GREY_FILL, edge=EDGE, lw=1.4, fs=9.5, bold=False, weight_lines=None, ha="center"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.12", fc=fill, ec=edge, lw=lw)
    ax.add_patch(p)
    ax.text(x + w / 2 if ha == "center" else x + 0.12, y + h / 2, text, ha=ha, va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.25)
    return p


def arrow(ax, p0, p1, color=EDGE, lw=1.3, ls="-", style="-|>", ms=11, conn="arc3,rad=0"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, linestyle=ls,
                                 connectionstyle=conn, shrinkA=0, shrinkB=0))


def poly(ax, pts, color=EDGE, lw=1.3, ls="-", head=True):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    if head:
        ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, ls=ls, solid_capstyle="butt")
        arrow(ax, pts[-2], pts[-1], color=color, lw=lw, ls=ls)
    else:
        ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="butt")


def save_tif(fig, name, dpi=1000):
    png = OUT / f"{name}.png"
    fig.savefig(png, dpi=dpi, facecolor="white")
    plt.close(fig)
    im = Image.open(png).convert("RGB")
    im.save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(dpi, dpi))
    return im.size


def figure1():
    fig = plt.figure(figsize=(7.48, 4.3))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 20.4)
    ax.set_ylim(0, 11.8)
    ax.axis("off")
    green, green_fill = "#2D6A57", "#E8F3EE"
    ax.text(0.7, 11.1, "A  Urban feedback", fontsize=10.5, fontweight="bold")
    ax.text(0.7, 10.4, "Nonlinear, cross-scale, path-dependent", fontsize=8.6, color=EDGE)
    ax.text(10.8, 11.1, "B  Plural learning", fontsize=10.5, fontweight="bold")
    ax.plot([10.15, 10.15], [5.0, 10.65], color="#D7DDE2", lw=0.8, zorder=0)

    box(ax, 3.1, 8.7, 3.8, 1.0, "New railway", fill=BLUE_FILL, edge=BLUE, fs=9.4, bold=True)
    box(ax, 0.8, 6.8, 3.9, 1.05, "Residents", fs=9.2)
    box(ax, 5.4, 6.8, 3.9, 1.05, "Land values", fs=9.2)
    box(ax, 0.8, 4.9, 3.9, 1.05, "Water demand", fs=9.2)
    box(ax, 5.4, 4.9, 3.9, 1.05, "Housing pressure", fs=9.2)
    box(ax, 2.0, 2.75, 6.1, 1.05, "New planning needs",
        fill=green_fill, edge=green, fs=9.3, bold=True)
    for start, end in (
        ((3.9, 8.7), (2.75, 7.85)), ((6.1, 8.7), (7.35, 7.85)),
        ((2.75, 6.8), (2.75, 5.95)), ((7.35, 6.8), (7.35, 5.95)),
        ((2.75, 4.9), (3.8, 3.8)), ((7.35, 4.9), (6.3, 3.8)),
    ):
        arrow(ax, start, end)

    for x, text in ((10.9, "Technical\nevidence"), (14.0, "Community\nexperience"),
                    (17.1, "Indigenous\nknowledge")):
        box(ax, x, 8.45, 2.8, 1.45, text, fs=8.8)
    box(ax, 11.2, 5.1, 8.3, 1.65, "Interpret together\nContext, coverage and validity",
        fill=BLUE_FILL, edge=BLUE, fs=9.1)
    for start, end in (
        ((12.3, 8.45), (13.0, 6.75)), ((15.4, 8.45), (15.35, 6.75)),
        ((18.5, 8.45), (17.7, 6.75)),
    ):
        arrow(ax, start, end)
    arrow(ax, (8.1, 3.275), (11.2, 5.9))
    box(ax, 12.4, 2.75, 5.6, 1.05, "Revised decisions",
        fill=green_fill, edge=green, fs=9.3, bold=True)
    arrow(ax, (15.35, 5.1), (15.2, 3.8), color=green)
    poly(ax, [(15.2, 2.75), (15.2, 1.55), (0.33, 1.55), (0.33, 9.2), (3.1, 9.2)],
         color=BLUE, ls=(0, (4, 2.5)))
    ax.text(8.1, 0.88, "Delayed feedback", fontsize=8.6, color=BLUE, ha="center")
    ax.text(10.2, 0.25, "Illustrative connections, not predictions",
            fontsize=8.4, color=EDGE, ha="center")
    return save_tif(fig, "Figure_1")


def figure2():
    fig = plt.figure(figsize=(7.48, 4.1))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 20.4)
    ax.set_ylim(0, 11.2)
    ax.axis("off")
    ax.text(0.3, 10.5, "Idealized approval-centred emphasis", fontsize=10, fontweight="bold")
    ax.text(0.3, 9.85, "Follow-up obligations vary by instrument and jurisdiction.", fontsize=8)
    xs = [0.3, 5.4, 10.5, 15.6]
    top = ["Criteria and\noptions", "Forecast and\nconsultation", "Decision and\nconditions",
           "Project delivery\nand follow-up"]
    bottom = ["Implementation\nand monitoring", "Accountable\nreview and action",
              "Recorded lessons\nand retained options", "Reassessment\nand next decisions"]
    for x, text in zip(xs, top):
        box(ax, x, 7.9, 4.5, 1.5, text, fs=9)
    for left, right in zip(xs, xs[1:]):
        arrow(ax, (left + 4.5, 8.65), (right, 8.65))
    ax.text(0.3, 6.85, "Conditional learning-oriented cycle", fontsize=10, fontweight="bold", color=BLUE)
    for x, text in zip(xs, bottom):
        box(ax, x, 4.35, 4.5, 1.8, text, fill=BLUE_FILL, edge=BLUE, fs=8.8)
    for left, right in zip(xs, xs[1:]):
        arrow(ax, (left + 4.5, 5.25), (right, 5.25), color=BLUE)
    poly(ax, [(17.85, 4.35), (17.85, 3.65), (2.55, 3.65), (2.55, 4.35)], color=BLUE)
    box(ax, 0.8, 0.35, 18.8, 2.45,
        "Independent capability and safeguard test\n"
        "Comparable baseline/follow-up + credible no-disturbance contrast + retention\n"
        "Performance and distributional safeguards over a declared horizon\n"
        "Design intent or operational correction alone does not verify capability gain.",
        fill="white", lw=1.8, fs=8.7)
    arrow(ax, (7.65, 4.35), (7.65, 2.8), ls=(0, (3, 2)))
    return save_tif(fig, "Figure_2")


# ---------------------------------------------------------------- Figure 3
def figure3():
    fig = plt.figure(figsize=(7.48, 5.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 20.4)
    ax.set_ylim(0, 14.2)
    ax.axis("off")
    fs, fi = 8.6, 7.6
    dash = (0, (4, 2.5))
    # row 2: disturbance chain, learning, response, performance
    y2, h2 = 8.3, 1.45
    box(ax, 0.15, y2, 1.95, h2, "Disturbance\n$\\mathbf{d}$", fs=fs)
    box(ax, 2.5, y2, 2.1, h2, "Exposure\n$\\mathbf{x}=\\mathbf{E}\\mathbf{d}$", fs=fs)
    box(ax, 5.0, y2, 2.6, h2, "Evidence\n$\\mathbf{e}=\\mathbf{M}\\mathbf{x}+\\mathbf{r}+\\boldsymbol{\\xi}$", fs=7.8)
    box(ax, 8.3, y2, 4.0, h2, "Decision-effective learning\n$\\ell=MQ\\cdot\\gamma\\cdot\\sigma$", fill=BLUE_FILL, edge=BLUE,
        lw=1.8, fs=7.9, bold=True)
    box(ax, 12.85, y2, 3.6, h2, "Response $\\mathbf{u}$\ndirection $(\\mathbf{A}_I+T\\mathbf{A}_T)\\hat{\\mathbf{e}}$\n$\\|\\mathbf{u}\\|\\leq\\kappa(F+O+R)$",
        fs=7.2)
    box(ax, 17.0, y2, 2.6, h2, "Sustainability\nperformance $\\mathbf{S}$", fill="white", lw=2.6, fs=fs + 0.4, bold=True)
    ym = y2 + h2 / 2
    for a, b in [(2.1, 2.5), (4.6, 5.0), (7.6, 8.3), (12.3, 12.85), (16.45, 17.0)]:
        arrow(ax, (a, ym), (b, ym))
    # row 1: learning inputs
    y1, h1 = 11.15, 1.5
    box(ax, 0.15, y1, 3.8, h1, "Valid routine input $\\mathbf{r}$\nretained in both arms", fill="white", edge=BLUE, fs=fi)
    arrow(ax, (3.95, 11.4), (5.05, y2 + h2), color=BLUE)
    box(ax, 4.5, y1, 4.05, h1, "Salience, attention threshold\n$\\sigma=\\|\\mathbf{e}\\|^h/(\\|\\mathbf{e}\\|^h+\\sigma_0^h)$",
        fill="white", edge=BLUE, fs=fi)
    box(ax, 8.85, y1, 3.75, h1, "Monitoring quality $MQ$\n(required-route veto;\nplurality scaled by $D$)", fill="white", edge=BLUE, fs=fi)
    box(ax, 12.9, y1, 2.4, h1, "Decision\nconnectivity $\\gamma$", fill="white", edge=BLUE, fs=fi)
    arrow(ax, (6.3, y2 + h2), (6.3, y1), color=BLUE)
    arrow(ax, (7.2, y1), (9.3, y2 + h2), color=BLUE)
    arrow(ax, (10.7, y1), (10.5, y2 + h2), color=BLUE)
    arrow(ax, (13.6, y1), (11.8, y2 + h2), color=BLUE)
    # losses and external change
    poly(ax, [(3.55, y2 + h2), (3.55, 10.1), (0.05, 10.1), (0.05, 13.55), (17.9, 13.55), (17.9, y2 + h2)])
    ax.text(10.7, 13.85, "losses from exposure ($-\\mathbf{B}\\mathbf{x}$)", ha="center", va="center", fontsize=fs)
    box(ax, 18.2, 11.35, 2.05, 1.15, "$\\boldsymbol{\\Omega}$: external\nchange", fill="white", edge="#888888", lw=1.0, fs=fi)
    arrow(ax, (19.2, 11.35), (19.2, y2 + h2))
    # row 3: validity, retention, retained learning, capacity update, capacity
    y3, h3 = 4.75, 1.55
    box(ax, 0.15, y3, 4.15, h3, "Latent validity of $\\mathbf{M}\\mathbf{x}+\\mathbf{r}$\nversus error $\\boldsymbol{\\xi}$;\naccuracy discards invalid lessons",
        fill="white", edge=BLUE, fs=7.2)
    box(ax, 5.0, y3, 2.6, h3, "Retention\n$\\tilde{\\rho}=\\rho\\,(\\theta+(1-\\theta)L)$", fill="white", edge=BLUE, fs=fi)
    box(ax, 8.3, y3, 4.0, h3, "Retained learning\n$\\ell^r=\\tilde{\\rho}\\,\\ell$", fill=BLUE_FILL, edge=BLUE, lw=1.8, fs=fs, bold=True)
    box(ax, 12.85, y3 - 0.1, 4.6, h3 + 0.2,
        "Capacity update, Eq. (2)\n$+\\,\\eta\\ell^r\\mathbf{q}\\odot$\n$[v(\\mathbf{1}-\\mathbf{C})-(1-v)(1-\\mu^{acc})\\mathbf{C}]$\n$-(\\delta_t\\mathbf{1}+\\boldsymbol{\\zeta}_t)\\odot\\mathbf{C}+\\boldsymbol{\\Omega}^C$",
        fs=7.0)
    box(ax, 17.75, y3, 2.2, h3, "Capacity $\\mathbf{C}$\n(proxy for\ncapability)", fill="white", lw=2.6, fs=fs, bold=True)
    ym3 = y3 + h3 / 2
    arrow(ax, (10.3, y2), (10.3, y3 + h3), color=BLUE)
    arrow(ax, (7.6, ym3), (8.3, ym3), color=BLUE)
    arrow(ax, (12.3, ym3), (12.85, ym3))
    arrow(ax, (17.45, ym3), (17.75, ym3))
    box(ax, 12.85, 6.55, 4.6, 1.15, "decay $\\delta_t=\\delta_0(1-\\varphi\\tilde{\\rho})$;\noption exercise and lock-in $\\boldsymbol{\\zeta}_t$",
        fill="white", edge="#888888", lw=1.0, fs=7.0)
    arrow(ax, (15.15, 6.55), (15.15, y3 + h3 + 0.1))
    poly(ax, [(5.6, y2), (5.6, 7.75), (2.2, 7.75), (2.2, y3 + h3)], color=BLUE)
    poly(ax, [(2.2, y3), (2.2, 4.05), (14.0, 4.05), (14.0, y3 - 0.1)], color=BLUE)
    # next-period feedback (dashed)
    poly(ax, [(18.4, y3 + h3), (18.4, 7.85), (14.65, 7.85), (14.65, y2)], ls=dash)
    poly(ax, [(18.7, y3), (18.7, 3.45), (6.3, 3.45), (6.3, y3)], ls=dash)
    poly(ax, [(4.75, 3.45), (4.75, 10.6), (10.0, 10.6), (10.0, y1)], ls=dash)
    ax.plot([4.75, 6.3], [3.45, 3.45], color=EDGE, lw=1.3, ls=dash)
    ax.text(12.4, 3.15, "dashed: acts in the next period", fontsize=7.2, ha="center", va="center", style="italic")
    # row 4: attribution, safeguards, classification
    y4, h4 = 0.25, 2.1
    box(ax, 0.15, y4, 7.9, h4, "Counterfactual period without disturbance\n($\\mathbf{d}=\\mathbf{0}$; same valid routine input):\nattributable gain $\\Delta^dC=\\boldsymbol{\\beta}^{\\top}(\\mathbf{C}_{t+1}-\\mathbf{C}^{0}_{t+1})$",
        fill="white", fs=fi)
    box(ax, 8.55, y4, 6.25, h4, "Safeguards: critical and temporary minima,\ngroup limits and cumulative budgets,\ncritical-capacity floors; no irreversible harm;\nmonitoring sufficiency $MQ\\geq MQ_{min}$",
        fill="white", fs=fi)
    box(ax, 15.3, y4, 4.3, h4, "Latent model class, Eq. (3)\nand capacity status\nEmpirical evidence:\nnot assessed in simulations",
        fs=7.3, bold=True)
    arrow(ax, (8.05, y4 + h4 / 2), (8.55, y4 + h4 / 2))
    arrow(ax, (14.8, y4 + h4 / 2), (15.3, y4 + h4 / 2))
    arrow(ax, (19.4, y3), (19.4, y4 + h4))
    poly(ax, [(19.6, ym), (20.25, ym), (20.25, y4 + h4 / 2), (19.6, y4 + h4 / 2)])
    return save_tif(fig, "Figure_3")

# ---------------------------------------------------------------- Figure 4
def figure4():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.48, 3.35), gridspec_kw=dict(width_ratios=[1, 1.25], wspace=0.32))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.16, top=0.80)
    fig.suptitle("Illustrative, uncalibrated model; not empirical validation", fontsize=9.5,
                 fontweight="bold", y=0.97)
    C0 = [0.4] * K
    S0 = [0.55, 0.60]
    xs = [i / 400 * 1.2 for i in range(401)]
    for h, col, ls, lab in [(2.0, "#1A1A1A", "-", "attention threshold ($h=2$)"), (1.0, BLUE, "--", "no threshold ($h=1$)")]:
        P = dict(BASE, h=h, noise=0.0)
        ys = []
        for dd in xs:
            P["sig"] = dd
            R = period(DESIGNS["B"], P, C0, S0, 1.0)
            ys.append(100 * (sum(R["Cn"]) - sum(C0)) / K)
        a1.plot(xs, ys, color=col, ls=ls, lw=1.8, label=lab)
    ystar = 0.4 * math.sqrt(1 / 3)
    a1.axvspan(0, ystar, color="#E3E3E3", lw=0)
    a1.text(ystar / 2, 1.85, "convex\ndomain", ha="center", va="top", fontsize=8, style="italic")
    a1.axhline(0, color="#999999", lw=0.8)
    a1.set_xlim(0, 1.2)
    a1.set_xlabel("Disturbance magnitude $\\|\\mathbf{d}\\|$", fontsize=9)
    a1.set_ylabel("Capacity change in period (×100)", fontsize=9)
    a1.legend(fontsize=7.5, frameon=False, loc="lower right")
    a1.set_title("(a) Fixed-state, fixed-direction response", fontsize=8.5, loc="left", fontweight="bold")
    a1.tick_params(labelsize=8)
    # panel b
    P = dict(BASE)
    ep = lambda t: t % 5 == 2
    hb = run_design("B", P, 30, ep)
    ha = run_design("A", P, 30, ep)
    long = run_design("B", P, 4000, ep)
    hm = run_design("B", P, 30, ep, C0=long["C"])
    t = list(range(31))
    for k in range(30):
        if ep(k):
            a2.axvline(k + 1, color="#E3E3E3", lw=2.2, zorder=0)
    a2.plot(t, hb["idx"], color="#1A1A1A", lw=1.8, label="Design B")
    a2.plot(t, hb["idxcf"], color="#1A1A1A", lw=1.2, ls=":", label="Design B, no disturbance")
    a2.plot(t, hm["idx"], color="#D55E00", lw=1.6, ls="-.", label="Design B, mature start")
    a2.plot(t, hm["idxcf"], color="#D55E00", lw=1.0, ls=":", label="mature start, no disturbance")
    a2.plot(t, ha["idx"], color=BLUE, lw=1.6, ls="--", label="Design A")
    a2.set_xlim(0, 30)
    a2.set_ylim(0, 0.5)
    a2.set_xlabel("Review period $t$ (years)", fontsize=9)
    a2.set_ylabel("Capacity index $C_t$", fontsize=9)
    a2.set_title("(b) Episodic latent capacity trajectories", fontsize=8.5, loc="left", fontweight="bold")
    a2.legend(fontsize=7, frameon=False, loc="lower left", ncol=1)
    a2.tick_params(labelsize=8)
    for a in (a1, a2):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    return save_tif(fig, "Figure_4")


# ---------------------------------------------------------------- Graphical abstract
def graphical_abstract():
    fig = plt.figure(figsize=(2657 / 300, 1063 / 300))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 26.57)
    ax.set_ylim(0, 10.63)
    ax.axis("off")
    ax.text(13.3, 9.95, "Disturbance-attributable capability gain", ha="center", va="center", fontsize=17,
            fontweight="bold")
    items = [("Experience", "Bounded shocks\nand observed\nimplementation", GREY_FILL, EDGE, 1.6),
             ("Learning", "Valid evidence\nreaches actors\nwith authority", BLUE_FILL, BLUE, 1.6),
             ("Retention", "Lessons recorded\nand reused;\noptions retained", BLUE_FILL, BLUE, 1.6),
             ("Capability test", "Routine-aware\ncontrast; retention\nand safeguards", "white", EDGE, 3.2)]
    x = 0.3
    for i, (title, body, fill, edge, lw) in enumerate(items):
        p = FancyBboxPatch((x, 4.55), 5.6, 4.35, boxstyle="round,pad=0,rounding_size=0.25", fc=fill, ec=edge, lw=lw)
        ax.add_patch(p)
        ax.text(x + 2.8, 8.25, title, ha="center", va="center", fontsize=16, fontweight="bold")
        ax.text(x + 2.8, 6.25, body, ha="center", va="center", fontsize=14.5, linespacing=1.3)
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + 5.72, 6.7), (x + 6.68, 6.7),
                                         arrowstyle="simple,head_width=1.2,head_length=1.0,tail_width=0.4",
                                         mutation_scale=20, fc="#1A1A1A", ec="none"))
        x += 6.75
    p = FancyBboxPatch((0.3, 0.3), 25.97, 3.7, boxstyle="round,pad=0,rounding_size=0.25", fc=GREY_FILL, ec=EDGE, lw=1.6)
    ax.add_patch(p)
    ax.text(13.3, 2.95, "Design is not verified gain: test attribution, retention and safeguards",
            ha="center", va="center", fontsize=14.5, fontweight="bold")
    ax.text(13.3, 2.05, "Convexity requires fixed state and direction within a bounded domain.",
            ha="center", va="center", fontsize=14.5)
    ax.text(13.3, 1.15, "Five cases inform mechanisms, not model validation. AI-assisted diagram code.",
            ha="center", va="center", fontsize=13.5)
    return save_tif(fig, "Graphical_Abstract", dpi=300)


if __name__ == "__main__":
    OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    OUT.mkdir(parents=True, exist_ok=True)
    check_against_script()
    export_model_records()
    print("Figure 1", figure1())
    print("Figure 2", figure2())
    print("Figure 3", figure3())
    print("Figure 4", figure4())
    print("Graphical abstract", graphical_abstract())
