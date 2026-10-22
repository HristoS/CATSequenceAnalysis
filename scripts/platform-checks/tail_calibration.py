"""Hold-out calibration of two tail models on the full score_null cells (psql output on stdin).

Each cell is split at random in two halves; a model is fitted on one and asked for the score whose
predicted p is 1e-3 (and 1e-2); the other half's exceedances of that score are counted. A calibrated
model expects 0.001 (0.01) of the held-out half above it.
"""
import math
import random
import sys

EULER = 0.5772156649015329


def gumbel_quantile(sample, p):
    mean = sum(sample) / len(sample)
    var = sum((x - mean) ** 2 for x in sample) / (len(sample) - 1)
    beta = math.sqrt(6 * var) / math.pi
    mu = mean - EULER * beta
    return mu - beta * math.log(-math.log(1 - p))


def exponential_tail_quantile(sample, p, tail=0.1):
    ordered = sorted(sample)
    u = ordered[int((1 - tail) * len(ordered))]
    excess = [x - u for x in ordered if x > u]
    sigma = sum(excess) / len(excess)
    # P(X > s) = tail * exp(-(s - u) / sigma)
    return u + sigma * math.log(tail / p)


def gpd_quantile(sample, p, tail=0.1):
    # Generalized Pareto over the top decile, method of moments: shape xi = (1 - m^2/v) / 2,
    # scale sigma = m (m^2/v + 1) / 2, with m and v the mean and variance of the excesses.
    ordered = sorted(sample)
    u = ordered[int((1 - tail) * len(ordered))]
    excess = [x - u for x in ordered if x > u]
    m = sum(excess) / len(excess)
    v = sum((e - m) ** 2 for e in excess) / (len(excess) - 1)
    xi = 0.5 * (1 - m * m / v)
    sigma = 0.5 * m * (m * m / v + 1)
    # P(X > s) = tail * (1 + xi (s - u) / sigma) ** (-1 / xi)
    return u + sigma / xi * ((p / tail) ** (-xi) - 1)


rng = random.Random(20261022)
totals = {}
cells = 0
for line in sys.stdin:
    parts = line.strip().split("|")
    if len(parts) != 3:
        continue
    source, window, raw = parts
    scores = [float(x) for x in raw.strip("{}").split(",")]
    if len(scores) < 10_000:
        continue
    cells += 1
    rng.shuffle(scores)
    fit, held = scores[:5_000], scores[5_000:]
    for name, quantile in (("gumbel", gumbel_quantile), ("exp-tail", exponential_tail_quantile), ("gpd", gpd_quantile)):
        for p in (1e-2, 1e-3, 1e-4):
            threshold = quantile(fit, p)
            key = (source, name, p)
            observed, expected = totals.get(key, (0, 0.0))
            totals[key] = (observed + sum(1 for x in held if x >= threshold), expected + p * len(held))

print(f"{cells} full cells")
print("source,model,p,observed,expected,ratio")
for (source, name, p), (observed, expected) in sorted(totals.items()):
    print(f"{source},{name},{p},{observed},{expected:.1f},{observed / expected:.2f}")
