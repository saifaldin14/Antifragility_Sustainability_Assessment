# Native verification and analysis script (Supplementary Material 1) for the article
# "Modelling Antifragility in Sustainability Assessment: A Protected Attribution Rule for Retained Capability Gain".
# Checks the internal consistency of the model and reproduces the illustrative, sensitivity, convexity, and
# horizon analyses reported in Supplementary Material 1. It does not validate the model against real
# project-governance systems.
# Run with: powershell -NoProfile -ExecutionPolicy Bypass -File run_verification.ps1
# Exit code 0: every test passed, every reported value was reproduced, and every test was evaluated on the
# reported number of cases; exit code 1 otherwise. Output does not depend on regional (number-format) settings.
param([switch]$InterfaceChecksOnly, [switch]$ScenarioChecksOnly, [switch]$DefinitionsOnly, [switch]$SafeguardChecksOnly, [switch]$BoundaryChecksOnly)
$ErrorActionPreference = 'Stop'
$src = @'
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

public static class AfVerify
{
    // Capacity positions: F flexibility, D diversity, O optionality, L learning capacity, R reversibility, T transformative capacity
    public const int iF = 0, iD = 1, iO = 2, iL = 3, iR = 4, iT = 5, K = 6;
    const int J = 15, m = 3, n = 4, p = 4, r = 3, NG = 3;
    public const double EpsS = 0.01, EpsC = 0.0, EpsD = 0.05, Eps0 = 0.05, MQmin = 0.10;
    const double TolBound = 1e-12, TolZero = 1e-15, TolConv = 1e-6;

    static Random Rng;
    static List<string> order;
    static Dictionary<string, long[]> res;
    static Dictionary<string, string> kind;
    static StringBuilder Out = new StringBuilder();
    static readonly System.Globalization.CultureInfo Inv = System.Globalization.CultureInfo.InvariantCulture;
    public static List<string> Counts = new List<string>();
    public static List<string> Reported = new List<string>();
    public static long TestFailures = 0;

    // ---------- utilities ----------
    static double U(double a, double b) { return a + (b - a) * Rng.NextDouble(); }
    static double[] Vec(int k, double a, double b) { var v = new double[k]; for (int i = 0; i < k; i++) v[i] = U(a, b); return v; }
    static double[] Simplex(int k) { var v = new double[k]; double s = 0; for (int i = 0; i < k; i++) { v[i] = -Math.Log(1 - Rng.NextDouble()); s += v[i]; } for (int i = 0; i < k; i++) v[i] /= s; return v; }
    static double[,] Mat(int a, int b, double lo, double hi) { var x = new double[a, b]; for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) x[i, j] = U(lo, hi); return x; }
    public static double[] Mul(double[,] A, double[] v)
    {
        if (A == null || v == null || A.GetLength(1) != v.Length)
            throw new ArgumentException("Matrix/vector dimensions do not agree.");
        foreach (double x in A) Finite(x, "matrix");
        Vector(v, v.Length, "vector", false, false);
        return MulChecked(A, v);
    }
    static double[] MulChecked(double[,] A, double[] v)
    {
        var output = new double[A.GetLength(0)];
        for (int i = 0; i < output.Length; i++)
        {
            double sum = 0; for (int j = 0; j < v.Length; j++) sum += A[i, j] * v[j];
            Finite(sum, "matrix product"); output[i] = sum;
        }
        return output;
    }
    public static double Dot(double[] a, double[] b)
    {
        if (a == null || b == null || a.Length != b.Length)
            throw new ArgumentException("Dot-product vectors must have matching dimensions.");
        Vector(a, a.Length, "left vector", false, false);
        Vector(b, b.Length, "right vector", false, false);
        return DotChecked(a, b);
    }
    static double DotChecked(double[] a, double[] b)
    {
        double sum = 0; for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        Finite(sum, "dot product"); return sum;
    }
    public static double Norm(double[] v)
    {
        if (v == null) throw new ArgumentNullException("v");
        double scale = 0, sum = 1;
        foreach (double value in v)
        {
            Finite(value, "norm vector"); double x = Math.Abs(value);
            if (x == 0) continue;
            if (x > scale) { double ratio = scale / x; sum = 1 + sum * ratio * ratio; scale = x; }
            else { double ratio = x / scale; sum += ratio * ratio; }
        }
        double norm = scale == 0 ? 0 : scale * Math.Sqrt(sum);
        Finite(norm, "norm result"); return norm;
    }
    public static double Clip(double x) { Finite(x, "projection input"); return Math.Max(0, Math.Min(1, x)); }
    static string F4(double x) { return x.ToString("F4", Inv); }
    static string F3(double x) { return x.ToString("F3", Inv); }

    // ---------- model ----------
    public class Inp
    {
        public double[] S, w, tau; public bool[] crit;
        public double[] C, beta, q;
        public double eta, delta0, chi, lockR, theta, phi, h, s0, kappa;
        public double[] d; public double[,] E, M; public double[] xi, routine, capacityFloor;
        public bool[] temporaryCritical;
        public bool cumulativeRequired;
        public double[] cumulativePrior, burdenIncrement, commonBurden, cumulativeLimit;
        public string monitoringForm = "product", salienceForm = "hill";
        public double mcov, macc, mtime, mplur, gamma, rho;
        public double[,] AI, AT, B, G; public double[] Om, OmC;
        public double[] pi; public double[,] Hm; public double Hmax; public double[] hgMax; public double xIrr;
        public double EQ, EC;
        public bool[] designatedGroups, reportingRouteRequired, reportingRouteAvailable;
        public Inp Copy() { return (Inp)MemberwiseClone(); }
    }
    public class Res
    {
        public double v, MQ, sigma, ell, ellR, rhoT, nu, deltaT, ne, bound; public double[] e, u, zeta, Cn, Sn, Smin, x;
        public double SC, SCn, SCmin, dC; public bool SB;
        public double[] groupHarm, cumulativeHarm;
        public double[] rawPerformance, rawImpactProxy, rawCapacity;
        public bool requiredRoutesPass;
    }
    public static double Mod(double theta, double c) { Probability(theta); Probability(c); return theta + (1 - theta) * c; }
    static double LogOnePlus(double x)
    {
        double shifted = 1 + x;
        return shifted == 1 ? x : Math.Log(shifted) * (x / (shifted - 1));
    }
    static double ExpMinusOne(double x)
    {
        double exponential = Math.Exp(x);
        if (exponential == 0) return -1;
        if (exponential == 1) return x;
        return (exponential - 1) * (x / Math.Log(exponential));
    }
    public static double ReviewRate(double annual, double interval)
    {
        Probability(annual); Positive(interval, "review interval");
        if (annual == 0 || annual == 1 || interval == 1) return annual;
        return -ExpMinusOne(interval * LogOnePlus(-annual));
    }
    static double AccurateSum(IEnumerable<double> values)
    {
        var parts = new List<double>();
        foreach (double value in values)
        {
            Finite(value, "summand");
            double high = value; int kept = 0;
            for (int j = 0; j < parts.Count; j++)
            {
                double low = parts[j];
                if (Math.Abs(high) < Math.Abs(low)) { double swap = high; high = low; low = swap; }
                double combined = high + low; Finite(combined, "sum");
                double residual = low - (combined - high);
                if (residual != 0) parts[kept++] = residual;
                high = combined;
            }
            parts.RemoveRange(kept, parts.Count - kept);
            parts.Add(high);
        }
        int next = parts.Count - 1;
        if (next < 0) return 0;
        double total = parts[next--], error = 0;
        while (next >= 0)
        {
            double value = parts[next--], combined = total + value;
            error = value - (combined - total); total = combined;
            if (error != 0) break;
        }
        if (next >= 0 && ((error < 0 && parts[next] < 0) || (error > 0 && parts[next] > 0)))
        {
            double correction = 2 * error, candidate = total + correction;
            if (candidate - total == correction) total = candidate;
        }
        Finite(total, "sum"); return total;
    }
    public static double Salience(double ne, double s0, double h)
    {
        Nonnegative(ne, "evidence norm"); Positive(s0, "salience scale"); Finite(h, "Hill exponent");
        if (h < 1) throw new ArgumentOutOfRangeException("h", "The Hill exponent must be at least one.");
        if (ne == 0) return 0;
        double logRatio = ne >= 0.5 * s0 && ne <= 2 * s0
            ? LogOnePlus((ne - s0) / s0) : Math.Log(ne) - Math.Log(s0);
        double z = h * logRatio;
        if (z >= 0) return 1 / (1 + Math.Exp(-z));
        double ez = Math.Exp(z); return ez / (1 + ez);
    }
    public static double Validity(double[] sig, double[] xi)
    {
        if (sig == null || xi == null || sig.Length != xi.Length)
            throw new ArgumentException("Signal and error require the same canonical evidence space.");
        double a = Norm(sig), b = Norm(xi), scale = Math.Max(a, b);
        return scale == 0 ? 0 : (a / scale) / (a / scale + b / scale);
    }

    static void Finite(double value, string name)
    {
        if (double.IsNaN(value) || double.IsInfinity(value))
            throw new ArgumentOutOfRangeException(name, "A finite value is required.");
    }
    static void Nonnegative(double value, string name) { Finite(value, name); if (value < 0) throw new ArgumentOutOfRangeException(name); }
    static void Positive(double value, string name) { Finite(value, name); if (value <= 0) throw new ArgumentOutOfRangeException(name); }
    static void Vector(double[] values, int size, string name, bool unit, bool nonnegative)
    {
        if (values == null || values.Length != size) throw new ArgumentException("Incorrect vector dimensions.", name);
        foreach (double value in values)
        {
            Finite(value, name);
            if ((unit || nonnegative) && value < 0 || unit && value > 1) throw new ArgumentOutOfRangeException(name);
        }
    }
    static void SimplexVector(double[] values, int size, string name)
    {
        Vector(values, size, name, true, true);
        if (size == 0 || Math.Abs(values.Sum() - 1) > 1e-12) throw new ArgumentException("Weights must sum to one.", name);
    }
    static void Mask(bool[] values, int size, string name)
    {
        if (values == null || values.Length != size) throw new ArgumentException("Explicit policy flags require matching dimensions.", name);
    }
    static void Matrix(double[,] values, int rows, int columns, string name, bool nonnegative)
    {
        if (values == null || values.GetLength(0) != rows || values.GetLength(1) != columns)
            throw new ArgumentException("Incorrect matrix dimensions.", name);
        foreach (double value in values) { Finite(value, name); if (nonnegative && value < 0) throw new ArgumentOutOfRangeException(name); }
    }
    public static void ValidateInput(Inp I)
    {
        if (I == null || I.S == null || I.S.Length == 0 || I.d == null || I.d.Length == 0 ||
            I.xi == null || I.xi.Length == 0 || I.E == null || I.E.GetLength(0) == 0 ||
            I.AI == null || I.AI.GetLength(0) == 0 || I.pi == null)
            throw new ArgumentException("Declare nonempty state, disturbance, exposure, evidence and response spaces.");
        int j = I.S.Length, m = I.d.Length, n = I.E.GetLength(0), p = I.xi.Length, r = I.AI.GetLength(0), g = I.pi.Length;
        Vector(I.S, j, "performance", true, true); SimplexVector(I.w, j, "performance weights");
        Vector(I.tau, j, "performance floors", true, true); Mask(I.crit, j, "critical performance");
        Mask(I.temporaryCritical, j, "mandatory temporary protections");
        Vector(I.C, K, "capacity", true, true); Vector(I.q, K, "relevance", true, true);
        SimplexVector(I.beta, K, "capacity weights"); Vector(I.capacityFloor, K, "capacity floors", true, true);
        foreach (double value in new[] { I.eta, I.delta0, I.chi, I.lockR, I.theta, I.phi, I.mcov, I.macc, I.mtime, I.mplur, I.gamma, I.rho, I.EQ, I.EC }) Probability(value);
        if (I.eta + I.delta0 + I.chi + I.lockR > 1)
            throw new ArgumentOutOfRangeException("rates", "S29 must hold for the actual review interval after rate conversion.");
        Positive(I.s0, "salience scale"); Positive(I.kappa, "response scale"); Finite(I.h, "Hill exponent");
        if (I.h < 1) throw new ArgumentOutOfRangeException("h");
        Vector(I.d, m, "disturbance", false, true); Matrix(I.E, n, m, "exposure", true);
        Matrix(I.M, p, n, "monitoring map", false); Vector(I.xi, p, "error", false, false);
        Vector(I.routine, p, "routine signal", false, false);
        Matrix(I.AI, r, p, "incremental response", false); Matrix(I.AT, r, p, "transformative response", false);
        Matrix(I.B, j, n, "performance loss", true); Matrix(I.G, j, r, "response effects", false);
        Vector(I.Om, j, "external performance change", false, false);
        if (I.OmC != null) Vector(I.OmC, K, "external capacity change", false, false);
        Matrix(I.Hm, g, n, "group harm map", true); Vector(I.pi, g, "equity weights", false, true);
        Vector(I.hgMax, g, "group harm limits", false, true); Mask(I.designatedGroups, g, "designated groups");
        Mask(I.reportingRouteRequired, g, "required group reporting routes");
        Mask(I.reportingRouteAvailable, g, "available group reporting routes");
        Nonnegative(I.Hmax, "aggregate harm limit"); Nonnegative(I.xIrr, "irreversibility threshold");
        if (I.monitoringForm != "product" && I.monitoringForm != "geometric" && I.monitoringForm != "minimum")
            throw new ArgumentException("Unknown monitoring aggregation form.");
        if (I.salienceForm != "hill" && I.salienceForm != "linear_capped") throw new ArgumentException("Unknown salience form.");
        if (I.cumulativeRequired)
        {
            Vector(I.cumulativePrior, g, "prior cumulative burdens", false, true);
            Vector(I.commonBurden, g, "common integrated burdens", false, true);
            Vector(I.burdenIncrement, g, "disturbance integrated burdens", false, true);
            Vector(I.cumulativeLimit, g, "cumulative budgets", false, true);
        }
    }
    public static bool RequiredRoutesHold(bool[] required, bool[] available)
    {
        if (required == null || available == null || required.Length != available.Length)
            throw new ArgumentException("Every declared group reporting route requires an explicit availability flag.");
        for (int g = 0; g < required.Length; g++) if (required[g] && !available[g]) return false;
        return true;
    }
    public static bool InitialProtectionHolds(Inp I)
    {
        ValidateInput(I);
        bool safe = CriticalHold(I.S, I.crit, I.tau) && CriticalHold(I.S, I.temporaryCritical, I.tau) && FloorsHold(I.C, I.capacityFloor);
        if (I.cumulativeRequired) safe &= CumulativeSafe(I.cumulativePrior, new double[I.pi.Length], I.cumulativeLimit);
        return safe;
    }
    public static bool ReviewProtectionHolds(Res R, double initialSC, double horizonTolerance)
    {
        Probability(initialSC); Probability(horizonTolerance);
        ValidateClassification(R.SB, R.MQ, R.SC, R.SCn, R.SCmin, 0, 0, R.dC);
        return R.SB && R.MQ >= MQmin && R.SCn >= R.SC - EpsS && R.SCn >= initialSC - horizonTolerance;
    }
    public static bool FloorsHold(double[] state, double[] floors)
    {
        if (state == null || floors == null || state.Length != floors.Length)
            throw new ArgumentException("Protected states and floors must have matching dimensions.");
        bool safe = true;
        for (int k = 0; k < state.Length; k++)
        {
            Finite(state[k], "state"); Finite(floors[k], "floors");
            if (state[k] < 0 || state[k] > 1 || floors[k] < 0 || floors[k] > 1)
                throw new ArgumentOutOfRangeException("floors", "Normalized states and floors must be within [0,1].");
            if (state[k] < floors[k]) safe = false;
        }
        return safe;
    }
    public static bool CriticalHold(double[] state, bool[] critical, double[] floors)
    {
        if (state == null || critical == null || floors == null || state.Length != critical.Length || state.Length != floors.Length)
            throw new ArgumentException("Critical outcomes, masks and thresholds must agree.");
        Vector(floors, floors.Length, "critical thresholds", true, true);
        var effective = new double[floors.Length];
        for (int j = 0; j < floors.Length; j++) effective[j] = critical[j] ? floors[j] : 0;
        return FloorsHold(state, effective);
    }
    public static bool CumulativeSafe(double[] prior, double[] increment, double[] limits)
    {
        if (prior == null || increment == null || limits == null || prior.Length != increment.Length || prior.Length != limits.Length)
            throw new ArgumentException("Cumulative burdens require matching prior, increment and budget vectors.");
        bool safe = true;
        for (int g = 0; g < prior.Length; g++)
        {
            Finite(prior[g], "prior"); Finite(increment[g], "increment"); Finite(limits[g], "limits");
            if (prior[g] < 0 || increment[g] < 0 || limits[g] < 0)
                throw new ArgumentOutOfRangeException("limits", "Burdens and budgets must be nonnegative.");
            double total = prior[g] + increment[g]; Finite(total, "cumulative burden");
            if (total > limits[g]) safe = false;
        }
        return safe;
    }
    public static double MonitoringQuality(double coverage, double timeliness, double plurality, string form)
    {
        var components = new[] { coverage, timeliness, plurality };
        foreach (double component in components)
        {
            Finite(component, "monitoring");
            if (component < 0 || component > 1) throw new ArgumentOutOfRangeException("monitoring");
        }
        double product = coverage * timeliness * plurality;
        if (form == "product") return product;
        if (form == "geometric") return Math.Pow(product, 1.0 / 3);
        if (form == "minimum") return components.Min();
        throw new ArgumentException("Unknown monitoring aggregation form.", "form");
    }
    public static double[] CanonicalEvidence(double[] values, double[] references, int[] channels, double[] weights, string[] sourceIds)
    {
        if (values == null || references == null || channels == null || weights == null || sourceIds == null || values.Length == 0 ||
            values.Length != references.Length || values.Length != channels.Length || values.Length != weights.Length || values.Length != sourceIds.Length)
            throw new ArgumentException("Canonical evidence requires matching observations, references, channels, fixed source weights and source/version IDs.");
        if (channels.Any(x => x < 0)) throw new ArgumentOutOfRangeException("channels");
        if (channels.Max() >= values.Length) throw new ArgumentException("Canonical channel identifiers must be contiguous from zero.");
        int count = channels.Max() + 1;
        var totals = new double[count]; var fused = new double[count]; var scales = new double[count];
        var unique = new Dictionary<string, int>(); var indices = new List<int>();
        for (int i = 0; i < values.Length; i++)
        {
            Finite(values[i], "values"); Finite(references[i], "references"); Finite(weights[i], "weights");
            if (references[i] <= 0 || weights[i] < 0) throw new ArgumentOutOfRangeException("references");
            if (string.IsNullOrWhiteSpace(sourceIds[i])) throw new ArgumentException("Stable original-source/version IDs are required.");
            double normalized = values[i] / references[i]; Finite(normalized, "normalized evidence");
            string key = channels[i].ToString(Inv) + ":" + sourceIds[i];
            int prior;
            if (unique.TryGetValue(key, out prior))
            {
                if (normalized != values[prior] / references[prior] || weights[i] != weights[prior])
                    throw new ArgumentException("Conflicting versions or weights for one source/channel must be reconciled before fusion.");
                continue;
            }
            unique.Add(key, i); indices.Add(i); scales[channels[i]] = Math.Max(scales[channels[i]], weights[i]);
        }
        for (int j = 0; j < count; j++)
        {
            if (scales[j] <= 0) throw new ArgumentException("Every contiguous canonical channel needs positive source weight.");
        }
        for (int j = 0; j < count; j++)
        {
            int channel = j;
            var members = indices.Where(i => channels[i] == channel).ToArray();
            totals[j] = AccurateSum(members.Select(i => weights[i] / scales[channel]));
            fused[j] = AccurateSum(members.Select(i =>
                ((weights[i] / scales[channel]) / totals[channel]) * (values[i] / references[i])));
        }
        foreach (double value in fused) Finite(value, "fused evidence");
        return fused;
    }
    public class RiskBounds
    {
        public readonly double aggregate, joint;
        public readonly double[] group;
        public RiskBounds(double aggregateUpper, double[] groupUpper, double jointUpper)
        {
            aggregate = aggregateUpper; group = groupUpper; joint = jointUpper;
        }
    }
    static void Probability(double value)
    {
        Finite(value, "probability");
        if (value < 0 || value > 1) throw new ArgumentOutOfRangeException("probability");
    }
    static void RiskTolerance(double value)
    {
        Probability(value); if (value == 1) throw new ArgumentOutOfRangeException("risk tolerance", "Risk tolerances must lie in [0,1), not [0,1].");
    }
    public static RiskBounds DiscreteHarmRisk(double[][] harms, double[] probabilities, double[] equityWeights,
                                             double aggregateLimit, double[] groupLimits, bool[] designated)
    {
        if (harms == null || probabilities == null || harms.Length == 0 || harms.Length != probabilities.Length ||
            equityWeights == null || groupLimits == null || designated == null ||
            equityWeights.Length != groupLimits.Length || designated.Length != groupLimits.Length)
            throw new ArgumentException("Harm scenarios, probabilities and group policies must have matching dimensions.");
        foreach (double probability in probabilities) Probability(probability);
        double total = AccurateSum(probabilities);
        if (Math.Abs(total - 1) > 1e-12) throw new ArgumentException("Scenario probabilities must sum to one.");
        Finite(aggregateLimit, "aggregateLimit");
        if (aggregateLimit < 0) throw new ArgumentOutOfRangeException("aggregateLimit");
        var aggregateMass = new List<double>(); var jointMass = new List<double>();
        var groupMass = new List<double>[groupLimits.Length];
        var groupRisk = new double[groupLimits.Length];
        for (int g = 0; g < groupLimits.Length; g++)
        {
            groupMass[g] = new List<double>();
            Finite(equityWeights[g], "equityWeights"); Finite(groupLimits[g], "groupLimits");
            if (equityWeights[g] < 0 || groupLimits[g] < 0) throw new ArgumentOutOfRangeException("groupLimits");
        }
        for (int s = 0; s < harms.Length; s++)
        {
            if (harms[s] == null || harms[s].Length != groupLimits.Length) throw new ArgumentException("Every harm scenario requires every group.");
            foreach (double harm in harms[s]) { Finite(harm, "harms"); if (harm < 0) throw new ArgumentOutOfRangeException("harms"); }
            double aggregateHarm = Dot(equityWeights, harms[s]); Finite(aggregateHarm, "aggregate harm");
            bool aggregateBreach = aggregateHarm > aggregateLimit, any = aggregateBreach;
            if (aggregateBreach) aggregateMass.Add(probabilities[s]);
            for (int g = 0; g < groupLimits.Length; g++)
                if (designated[g] && harms[s][g] > groupLimits[g]) { groupMass[g].Add(probabilities[s]); any = true; }
            if (any) jointMass.Add(probabilities[s]);
        }
        for (int g = 0; g < groupRisk.Length; g++) groupRisk[g] = AccurateSum(groupMass[g]) / total;
        return new RiskBounds(AccurateSum(aggregateMass) / total, groupRisk, AccurateSum(jointMass) / total);
    }
    public static bool RiskAccepts(RiskBounds risk, double aggregateMaximum, double[] groupMaximum, double jointMaximum)
    {
        if (risk == null || risk.group == null || groupMaximum == null || risk.group.Length != groupMaximum.Length)
            throw new ArgumentException("Every designated group needs a supplied risk bound and risk limit.");
        Probability(risk.aggregate); Probability(risk.joint); RiskTolerance(aggregateMaximum); RiskTolerance(jointMaximum);
        bool safe = risk.aggregate <= aggregateMaximum && risk.joint <= jointMaximum;
        for (int g = 0; g < groupMaximum.Length; g++)
        {
            Probability(risk.group[g]); RiskTolerance(groupMaximum[g]);
            if (risk.group[g] > groupMaximum[g]) safe = false;
        }
        return safe;
    }
    public static bool FeasibleAlternative(double[] performance, bool[] critical, double[] thresholds,
                                          double[] capacity, double[] capacityFloors, bool irreversible, bool mandatoryPathSafe,
                                          RiskBounds risk, double aggregateMaximum, double[] groupMaximum, double jointMaximum)
    {
        bool performanceSafe = CriticalHold(performance, critical, thresholds);
        bool capacitySafe = FloorsHold(capacity, capacityFloors);
        bool riskSafe = RiskAccepts(risk, aggregateMaximum, groupMaximum, jointMaximum);
        return performanceSafe && capacitySafe && !irreversible && mandatoryPathSafe && riskSafe;
    }
    // Capacity update, shared by every test and analysis
    public static double[] UpdateC(double[] C, double eta, double ellR, double[] q, double v, double acc, double deltaT, double[] zeta, double[] OmC)
    {
        if (C == null || C.Length == 0) throw new ArgumentException("A nonempty capacity vector is required.");
        Vector(C, C.Length, "capacity", true, true); Vector(q, C.Length, "relevance", true, true);
        Vector(zeta, C.Length, "consumption rates", true, true);
        foreach (double value in new[] { eta, ellR, v, acc, deltaT }) Probability(value);
        if (eta + deltaT + zeta.Max() > 1) throw new ArgumentOutOfRangeException("rates", "The sufficient componentwise capacity-domain bound must hold.");
        if (OmC != null) Vector(OmC, C.Length, "external capacity change", false, false);
        return UpdateCChecked(C, eta, ellR, q, v, acc, deltaT, zeta, OmC);
    }
    static double[] UpdateCChecked(double[] C, double eta, double ellR, double[] q, double v, double acc, double deltaT, double[] zeta, double[] OmC, double[] raw = null)
    {
        var o = new double[C.Length];
        for (int k = 0; k < C.Length; k++)
        {
            double learning = eta * ellR * q[k];
            double invalidation = learning * (1 - v) * (1 - acc);
            double budget = eta + deltaT + zeta[k];
            if (budget > 1 + 4 * 2.2204460492503131e-16)
                throw new ArithmeticException("The capacity survival coefficient violates the admitted rate bound.");
            // S29 bounds the coefficient, not the stock; absorb only coefficient-summation roundoff.
            double survival = (1 - Math.Min(1, budget)) + (eta - invalidation);
            double z = learning * v * (1 - C[k]) + survival * C[k] + (OmC == null ? 0 : OmC[k]);
            Finite(z, "raw capacity update"); if (raw != null) raw[k] = z; o[k] = OmC == null ? z : Clip(z);
        }
        return o;
    }

    // One assessment period; zeroD evaluates the no-disturbance counterfactual on the same inputs
    public static Res Period(Inp I, bool zeroD)
    {
        ValidateInput(I);
        return PeriodChecked(I, zeroD);
    }
    static Res PeriodChecked(Inp I, bool zeroD)
    {
        var R = new Res(); var C = I.C;
        var d = zeroD ? new double[I.d.Length] : I.d;
        R.x = MulChecked(I.E, d); var sig = MulChecked(I.M, R.x);
        if (I.routine == null || I.routine.Length != sig.Length)
            throw new ArgumentException("Declare a valid routine signal for every canonical evidence channel.");
        for (int i = 0; i < sig.Length; i++) { Finite(I.routine[i], "routine"); sig[i] += I.routine[i]; }
        R.e = new double[sig.Length]; for (int i = 0; i < sig.Length; i++) R.e[i] = sig[i] + I.xi[i];
        R.v = Validity(sig, I.xi); R.ne = Norm(R.e);
        R.MQ = MonitoringQuality(I.mcov, I.mtime, I.mplur * Mod(I.theta, C[iD]), I.monitoringForm);
        R.requiredRoutesPass = RequiredRoutesHold(I.reportingRouteRequired, I.reportingRouteAvailable);
        if (!R.requiredRoutesPass) R.MQ = 0;
        if (I.salienceForm == "hill") R.sigma = Salience(R.ne, I.s0, I.h);
        else if (I.salienceForm == "linear_capped") R.sigma = Math.Min(1, R.ne / I.s0);
        else throw new ArgumentException("Unknown salience form.");
        R.ell = R.MQ * I.gamma * R.sigma;
        R.rhoT = I.rho * Mod(I.theta, C[iL]);
        R.ellR = R.rhoT * R.ell;
        int rr = I.AI.GetLength(0), pp = I.AI.GetLength(1);
        var eh = new double[pp]; if (R.ne > 0) for (int i = 0; i < pp; i++) eh[i] = R.e[i] / R.ne;
        var A = new double[rr, pp]; for (int i = 0; i < rr; i++) for (int j = 0; j < pp; j++) A[i, j] = I.AI[i, j] + C[iT] * I.AT[i, j];
        var us = MulChecked(A, eh); for (int i = 0; i < rr; i++) us[i] *= R.ell;
        R.bound = I.kappa * (C[iF] + C[iO] + C[iR]); Finite(R.bound, "response bound"); double nus = Norm(us);
        if (nus > R.bound) for (int i = 0; i < rr; i++) us[i] = R.bound > 0 ? (us[i] / nus) * R.bound : 0;
        R.u = us; R.nu = R.bound > 0 ? Math.Min(1, Norm(us) / R.bound) : 0;
        int JJ = I.S.Length; var Bx = MulChecked(I.B, R.x); var Gu = MulChecked(I.G, R.u);
        R.Smin = new double[JJ]; R.Sn = new double[JJ];
        R.rawImpactProxy = new double[JJ]; R.rawPerformance = new double[JJ];
        for (int j = 0; j < JJ; j++)
        {
            R.rawImpactProxy[j] = I.S[j] - Bx[j]; R.rawPerformance[j] = R.rawImpactProxy[j] + Gu[j] + I.Om[j];
            R.Smin[j] = Clip(R.rawImpactProxy[j]); R.Sn[j] = Clip(R.rawPerformance[j]);
        }
        R.SC = DotChecked(I.w, I.S); R.SCn = DotChecked(I.w, R.Sn); R.SCmin = DotChecked(I.w, R.Smin);
        R.deltaT = I.delta0 * (1 - I.phi * R.rhoT);
        R.zeta = new double[C.Length]; R.zeta[iO] = I.chi * R.nu; R.zeta[iR] = I.chi * R.nu + I.lockR;
        R.rawCapacity = new double[C.Length];
        R.Cn = UpdateCChecked(C, I.eta, R.ellR, I.q, R.v, I.macc, R.deltaT, R.zeta, I.OmC, R.rawCapacity);
        R.dC = DotChecked(I.beta, R.Cn) - DotChecked(I.beta, C);
        // Safeguards: critical thresholds, aggregate and designated group harm limits, irreversibility
        var hg = MulChecked(I.Hm, R.x); double H = DotChecked(I.pi, hg);
        bool sb = H <= I.Hmax;
        R.groupHarm = hg;
        for (int g = 0; g < hg.Length; g++) if (I.designatedGroups[g] && hg[g] > I.hgMax[g]) sb = false;
        if (R.x.Length > 0 && R.x.Max() > I.xIrr) sb = false;
        if (!CriticalHold(R.Sn, I.crit, I.tau) || !FloorsHold(R.Cn, I.capacityFloor)) sb = false;
        if (!CriticalHold(R.Smin, I.temporaryCritical, I.tau)) sb = false;
        if (I.cumulativeRequired)
        {
            if (I.commonBurden == null || I.burdenIncrement == null || I.commonBurden.Length != hg.Length || I.burdenIncrement.Length != hg.Length)
                throw new ArgumentException("Declare both common and disturbance-linked integrated burdens for every cumulative receptor.");
            for (int g = 0; g < hg.Length; g++)
            {
                Finite(I.commonBurden[g], "common burden"); Finite(I.burdenIncrement[g], "disturbance burden");
                if (I.commonBurden[g] < 0 || I.burdenIncrement[g] < 0)
                    throw new ArgumentOutOfRangeException("burden", "Integrated burdens must be nonnegative.");
            }
            var increment = I.commonBurden.Select((value, g) => value + (zeroD ? 0 : I.burdenIncrement[g])).ToArray();
            if (!CumulativeSafe(I.cumulativePrior, increment, I.cumulativeLimit)) sb = false;
            R.cumulativeHarm = I.cumulativePrior.Select((value, g) => value + increment[g]).ToArray();
        }
        R.SB = sb;
        ValidateClassification(R.SB, R.MQ, R.SC, R.SCn, R.SCmin, Norm(d), 0, R.dC);
        return R;
    }

    // These are latent model classes, not empirically verified case classifications.
    public static readonly string[] LatentClassNames = { "Unacceptable", "Unverified", "Fragile", "CapabilityGain", "Robust", "Resilient" };
    public static int LatentClassOf(bool SB, double MQ, double SC, double SCn, double SCmin, double nd, double dCd, double dC)
    {
        ValidateClassification(SB, MQ, SC, SCn, SCmin, nd, dCd, dC);
        if (!SB) return 0;
        if (MQ < MQmin) return 1;
        if (SCn < SC - EpsS) return 2;
        if (nd > EpsD && dCd > EpsC && dC >= -EpsC) return 3;
        return SCmin >= SC - EpsS ? 4 : 5;
    }
    static void ValidateClassification(bool SB, double MQ, double SC, double SCn, double SCmin, double nd, double dCd, double dC)
    {
        foreach (double value in new[] { MQ, SC, SCn, SCmin }) Probability(value);
        Nonnegative(nd, "disturbance norm");
        foreach (double value in new[] { dCd, dC })
        {
            Finite(value, "capacity contrast");
            if (value < -1 || value > 1) throw new ArgumentOutOfRangeException("capacity contrast");
        }
    }
    public const string SimulationEvidenceStatus = "NotAssessed_IllustrativeSimulation";
    public static string CapacityEvidenceStatus(double EC, bool isSimulation, bool[] established)
    {
        if (double.IsNaN(EC) || double.IsInfinity(EC) || EC < 0 || EC > 1)
            throw new ArgumentOutOfRangeException("EC", "Evidence quality must be finite and within [0,1].");
        if (established == null || established.Length != 8)
            throw new ArgumentException("Eight G1-G8 evidence criteria are required.", "established");
        if (isSimulation) return SimulationEvidenceStatus;
        if (EC == 0 || established.Any(x => !x)) return "InsufficientEvidence";
        return "ReadyForIndependentEffectAndSafeguardAnalysis";
    }
    public static int Status(double dC) { Finite(dC, "capacity change"); if (Math.Abs(dC) > 1) throw new ArgumentOutOfRangeException("capacity change"); return dC > EpsC ? 1 : (dC < -EpsC ? -1 : 0); }
    public static double Cgi(double dCd, double EC, double EQ, bool SB, double nd)
    {
        Finite(dCd, "capacity contrast"); if (Math.Abs(dCd) > 1) throw new ArgumentOutOfRangeException("capacity contrast");
        Probability(EC); Probability(EQ); Nonnegative(nd, "disturbance norm");
        double psi = dCd > 0 ? EC * EQ * (SB ? 1 : 0) : 1;
        return psi * dCd / (nd + Eps0);
    }
    public static int SelectBest(double[] V, bool[] feas)
    {
        if (V == null || feas == null || V.Length != feas.Length) throw new ArgumentException("Alternative values and feasibility flags must agree.");
        int best = -1;
        for (int i = 0; i < V.Length; i++) { Finite(V[i], "values"); if (feas[i] && (best < 0 || V[i] > V[best])) best = i; }
        return best;
    }

    // ---------- bookkeeping ----------
    static void Check(string name, string k, bool ok)
    {
        if (!res.ContainsKey(name)) { res[name] = new long[2]; order.Add(name); kind[name] = k; }
        res[name][0]++; if (!ok) res[name][1]++;
    }
    static void Reset() { order = new List<string>(); res = new Dictionary<string, long[]>(); kind = new Dictionary<string, string>(); }
    static bool IsTest(string k) { return k == "Property" || k == "Implementation" || k == "Consistency" || k == "Edge"; }
    static string Id(string name) { return name.Split(' ')[0]; }
    static void Dump(string header, string tag)
    {
        Out.AppendLine("== " + header);
        foreach (var nm in order)
        {
            long nObs = res[nm][0], fail = res[nm][1];
            Out.AppendLine(string.Format(Inv, "{0}|{1}|N={2}|violations={3}", kind[nm], nm, nObs, fail));
            if (IsTest(kind[nm])) TestFailures += fail;
            if (IsTest(kind[nm]) && kind[nm] != "Edge") Counts.Add(tag + "|" + Id(nm) + "|N=" + nObs);
            else if (!IsTest(kind[nm])) Reported.Add(tag + "|" + Id(nm) + "|N=" + nObs + "|count=" + fail);
        }
    }
    static void Rep(string key, string val) { Reported.Add(key + "|" + val); Out.AppendLine(key + "|" + val); }

    // ---------- randomized scenarios ----------
    class Regime { public string Name, Tag; public int Seed; public double SLo, TauHi, EtaLo, DeltaCap, MQLo, RGLo, BLo, BHi, GLo, GHi, HmLo, HmHi, hgLo, hgHi, irLo, irHi; }
    static readonly Regime RA = new Regime { Name = "A (stress)", Tag = "A", Seed = 20261006, SLo = 0, TauHi = 0.3, EtaLo = 0, DeltaCap = 1, MQLo = 0, RGLo = 0, BLo = 0.15, BHi = 0.15, GLo = -0.2, GHi = 0.5, HmLo = 0.5, HmHi = 6, hgLo = 0.2, hgHi = 2, irLo = 1, irHi = 5 };
    static readonly Regime RB = new Regime { Name = "B (coverage)", Tag = "B", Seed = 20261007, SLo = 0.2, TauHi = 0.15, EtaLo = 0.2, DeltaCap = 0.15, MQLo = 0.6, RGLo = 0.3, BLo = 0.001, BHi = 0.12, GLo = -0.1, GHi = 0.6, HmLo = 1, HmHi = 10, hgLo = 0.5, hgHi = 3, irLo = 2, irHi = 8 };

    static Inp Draw(Regime R)
    {
        var I = new Inp();
        I.S = Vec(J, R.SLo, 1); I.w = Simplex(J); I.tau = Vec(J, 0, R.TauHi);
        I.crit = new bool[J]; for (int j = 0; j < J; j++) I.crit[j] = Rng.NextDouble() < 0.3;
        I.C = Vec(K, 0, 1); I.beta = Simplex(K); I.q = Vec(K, 0, 1);
        I.eta = U(R.EtaLo, 1); I.delta0 = U(0, Math.Min(R.DeltaCap, 1 - I.eta));
        double rem = 1 - I.eta - I.delta0; I.chi = U(0, Math.Min(0.3, rem)); I.lockR = U(0, Math.Min(0.05, rem - I.chi));
        I.theta = Rng.NextDouble() < 0.2 ? 1 : U(0, 1); I.phi = U(0, 1); I.h = U(1, 3); I.s0 = U(0.1, 2); I.kappa = U(0, 1);
        I.d = Rng.NextDouble() < 0.05 ? new double[m] : Vec(m, 0, 1.5);
        I.E = Mat(n, m, 0, 1); I.M = Mat(p, n, 0, 1.5); I.xi = Vec(p, -0.2, 0.2);
        I.routine = new double[p]; I.capacityFloor = new double[K]; I.temporaryCritical = new bool[J];
        I.cumulativeRequired = false;
        I.mcov = U(R.MQLo, 1); I.macc = U(R.MQLo, 1); I.mtime = U(R.MQLo, 1); I.mplur = U(R.MQLo, 1);
        I.rho = Rng.NextDouble() < 0.1 ? 0 : U(R.RGLo, 1); I.gamma = Rng.NextDouble() < 0.1 ? 0 : U(R.RGLo, 1);
        I.AI = Mat(r, p, -1, 1); I.AT = Mat(r, p, -1, 1);
        double bHi = U(R.BLo, R.BHi); I.B = Mat(J, n, 0, bHi); I.G = Mat(J, r, R.GLo, R.GHi); I.Om = Vec(J, -0.05, 0.05);
        I.OmC = Rng.NextDouble() < 0.5 ? null : Vec(K, -0.05, 0.05);
        I.pi = Vec(NG, 0.5, 2); I.Hm = Mat(NG, n, 0, 0.5); I.Hmax = U(R.HmLo, R.HmHi);
        I.hgMax = new double[NG]; I.designatedGroups = new bool[NG];
        for (int g = 0; g < NG; g++) { I.designatedGroups[g] = Rng.NextDouble() < 0.5; I.hgMax[g] = I.designatedGroups[g] ? U(R.hgLo, R.hgHi) : 0; }
        I.reportingRouteRequired = Enumerable.Repeat(true, NG).ToArray();
        I.reportingRouteAvailable = Enumerable.Repeat(true, NG).ToArray();
        I.xIrr = U(R.irLo, R.irHi); I.EQ = U(0, 1); I.EC = U(0, 1);
        return I;
    }

    static void RunRegime(Regime R, int trials)
    {
        Rng = new Random(R.Seed); Reset();
        long[] cls = new long[6]; long[,] st = new long[6, 3];
        for (int t = 0; t < trials; t++)
        {
            var I = Draw(R);
            var A = Period(I, false); var Z = Period(I, true);
            double nd = Norm(I.d), dCd = Dot(I.beta, A.Cn) - Dot(I.beta, Z.Cn);

            Check("T1 Learning lies in [0,1] and retained learning does not exceed it", "Property", A.ell >= 0 && A.ell <= 1 && A.ellR >= 0 && A.ellR <= A.ell);
            Check("T2 Implemented response respects kappa(F+O+R)", "Implementation", Norm(A.u) <= A.bound + TolBound);
            Check("T3 Performance remains in [0,1]", "Implementation", A.Sn.All(s => s >= 0 && s <= 1));
            if (I.OmC == null)
            {
                bool cOk = true, tight = true;
                for (int k = 0; k < K; k++)
                {
                    double a = I.eta * A.ellR * I.q[k], lo = a * A.v, hi = 1 - A.deltaT - A.zeta[k] - a * (1 - A.v) * (1 - I.macc);
                    if (A.Cn[k] < -TolBound || A.Cn[k] > 1 + TolBound) cOk = false;
                    if (A.Cn[k] < lo - TolBound || A.Cn[k] > hi + TolBound) tight = false;
                }
                Check("T4 Capacity remains in [0,1] without truncation (Property 1)", "Property", cOk);
                Check("T5 Capacity lies within the tight bounds of Property 1", "Property", tight);
                var qv = new double[K]; var dz = new double[K];
                for (int k = 0; k < K; k++) { qv[k] = I.q[k] * (A.v * (1 - I.C[k]) - (1 - A.v) * (1 - I.macc) * I.C[k]); dz[k] = (A.deltaT + A.zeta[k]) * I.C[k]; }
                double lhs = I.eta * A.ellR * Dot(I.beta, qv), rhs = Dot(I.beta, dz);
                if (Math.Abs(lhs - rhs) > TolBound) Check("T6 Accounting identity predicts the sign of capacity change (Property 3)", "Property", (A.dC > 0) == (lhs > rhs));
                if (I.rho == 0 || I.gamma == 0) Check("T7 No capacity gain without retention or connectivity (Property 2)", "Property", A.ellR == 0 && A.dC <= TolZero);
                Check("T8 No capacity gain from evidence without signal (Property 2)", "Property", Z.v == 0 && Z.dC <= TolZero);
                double bump = U(0, 1 - A.ellR), dC2 = Dot(I.beta, UpdateC(I.C, I.eta, A.ellR + bump, I.q, A.v, I.macc, A.deltaT, A.zeta, null)) - Dot(I.beta, I.C);
                double slope = I.eta * Dot(I.beta, qv);
                Check("T10 More retained learning raises capacity only when its net validity term is positive (Property 9)", "Property", slope >= 0 ? dC2 >= A.dC - TolBound : dC2 <= A.dC + TolBound);
            }
            if (nd == 0) Check("T9 No disturbance-attributable gain without disturbance", "Property", Math.Abs(dCd) <= TolZero);
            // Property 9: learning is non-decreasing in its inputs, including diversity and learning capacity
            var I2 = I.Copy(); I2.mcov += U(0, 1 - I.mcov); I2.mplur += U(0, 1 - I.mplur); I2.gamma += U(0, 1 - I.gamma); I2.rho += U(0, 1 - I.rho);
            I2.C = (double[])I.C.Clone(); I2.C[iD] += U(0, 1 - I.C[iD]); I2.C[iL] += U(0, 1 - I.C[iL]);
            double sc = U(1, 3); var B2 = Period(I2, false);
            double ell2 = B2.MQ * I2.gamma * Salience(A.ne * sc, I.s0, I.h), ellR2 = B2.rhoT * ell2;
            Check("T11 Learning is non-decreasing in monitoring, connectivity, retention, evidence strength, diversity, and learning capacity (Property 9)", "Property", ell2 >= A.ell - TolZero && ellR2 >= A.ellR - TolZero);

            bool notDown = A.SCn >= A.SC - EpsS, af = nd > EpsD && dCd > EpsC && A.dC >= -EpsC;
            bool ok = A.SB && A.MQ >= MQmin;
            bool[] k6 = { !A.SB, A.SB && A.MQ < MQmin, ok && !notDown, ok && notDown && af, ok && notDown && !af && A.SCmin >= A.SC - EpsS, ok && notDown && !af && A.SCmin < A.SC - EpsS };
            int c = LatentClassOf(A.SB, A.MQ, A.SC, A.SCn, A.SCmin, nd, dCd, A.dC);
            Check("T12 Exactly one response class is assigned, and the classification function agrees (Property 6)", "Property", k6.Count(z => z) == 1 && k6[c]);
            cls[c]++; st[c, Status(A.dC) + 1]++;
            bool strong = A.SB && A.MQ >= MQmin && notDown && af;
            Check("T13 Strong condition coincides with the capability-gain class", "Consistency", strong == (c == 3));
            double cgi = Cgi(dCd, I.EC, I.EQ, A.SB, nd);
            Check("T14 Index sign is consistent with attributable capacity change (Property 7)", "Property", ((cgi > 0) == (dCd > 0 && A.SB && I.EC * I.EQ > 0)) && ((cgi < 0) == (dCd < 0)));
            Check("T15 Index is bounded by 1/eps0 (Property 10)", "Property", Math.Abs(cgi) <= 1 / Eps0 + TolBound);
            if (c == 3) Check("T16 Capability-gain class implies a positive index when EC*EQ > 0", "Consistency", I.EC * I.EQ <= 0 || cgi > 0);
        }
        Dump("REGIME " + R.Name + " (seed " + R.Seed + ", " + trials + " scenarios)", R.Tag);
        Rep(R.Tag + "|classes", string.Join(",", cls));
        var stl = new List<string>(); for (int i = 0; i < 6; i++) stl.Add(st[i, 2] + "/" + st[i, 1] + "/" + st[i, 0]);
        Rep(R.Tag + "|status(gaining/stable/eroding by class)", string.Join(",", stl));
    }

    // ---------- structural properties ----------
    static void RunStructure()
    {
        Rng = new Random(20261012); Reset();
        for (int t = 0; t < 20000; t++)
        {
            var C = Vec(K, 0, 1); var beta = Simplex(K); var q = Vec(K, 0, 1); var zeta = Vec(K, 0, 0.1);
            double eta = U(0, 0.8), dl = U(0, 0.1), ellR = U(0, 1), v = U(0, 1), acc = U(0, 1); int k = Rng.Next(K); double s = U(0, 1);
            var C2 = C.Concat(new[] { C[k] }).ToArray(); var b2 = beta.Concat(new[] { beta[k] * (1 - s) }).ToArray(); b2[k] = beta[k] * s;
            var q2 = q.Concat(new[] { q[k] }).ToArray(); var z2 = zeta.Concat(new[] { zeta[k] }).ToArray();
            double g1 = Dot(beta, UpdateC(C, eta, ellR, q, v, acc, dl, zeta, null)) - Dot(beta, C);
            double g2 = Dot(b2, UpdateC(C2, eta, ellR, q2, v, acc, dl, z2, null)) - Dot(b2, C2);
            Check("T17 Capacity index is invariant to splitting a capacity into identical parts (Property 11)", "Property", Math.Abs(g1 - g2) <= TolBound);
        }
        for (int t = 0; t < 20000; t++)
        {
            var I = Draw(RB); I.theta = 1; I.phi = 0; I.chi = 0; I.lockR = 0; I.h = 1; I.OmC = null; I.C[iF] = 0; I.AT = new double[r, p]; I.macc = 1; I.xi = new double[p];
            var A = Period(I, false);
            double ne = A.ne, ell = I.mcov * I.macc * I.mtime * I.mplur * I.gamma * ne / (ne + I.s0), ellR = I.rho * ell;   // reduced form: four-component MQ, unbounded-validity update
            var Cr = new double[K]; for (int k = 0; k < K; k++) Cr[k] = I.C[k] + I.eta * ellR * I.q[k] * (1 - I.C[k]) - I.delta0 * I.C[k];
            bool ok = Math.Abs(A.ell - ell) <= TolBound && Enumerable.Range(0, K).All(k => Math.Abs(A.Cn[k] - Cr[k]) <= TolBound) && Math.Abs(A.bound - I.kappa * (I.C[iO] + I.C[iR])) <= TolBound;
            Check("T18 Switching off the extensions recovers the reduced form", "Implementation", ok);
        }
        Dump("STRUCTURAL PROPERTIES (seed 20261012)", "P");
    }

    // ---------- convergence, feedback equilibria, capability trap ----------
    class FeedbackMap
    {
        readonly Inp I;
        readonly double[] signal, ai, at;
        readonly double validity, salience;
        public FeedbackMap(Inp input)
        {
            ValidateInput(input); I = input;
            signal = MulChecked(I.M, MulChecked(I.E, I.d));
            var evidence = new double[signal.Length];
            for (int j = 0; j < signal.Length; j++) { signal[j] += I.routine[j]; evidence[j] = signal[j] + I.xi[j]; }
            validity = Validity(signal, I.xi); double norm = Norm(evidence);
            salience = I.salienceForm == "hill" ? Salience(norm, I.s0, I.h) : Math.Min(1, norm / I.s0);
            var direction = new double[evidence.Length];
            if (norm > 0) for (int j = 0; j < direction.Length; j++) direction[j] = evidence[j] / norm;
            ai = MulChecked(I.AI, direction); at = MulChecked(I.AT, direction);
        }
        public double[] Step(double[] C)
        {
            double mq = MonitoringQuality(I.mcov, I.mtime, I.mplur * Mod(I.theta, C[iD]), I.monitoringForm);
            if (!RequiredRoutesHold(I.reportingRouteRequired, I.reportingRouteAvailable)) mq = 0;
            double ell = mq * I.gamma * salience, rho = I.rho * Mod(I.theta, C[iL]);
            var desired = new double[ai.Length];
            for (int j = 0; j < desired.Length; j++) desired[j] = ell * (ai[j] + C[iT] * at[j]);
            double bound = I.kappa * (C[iF] + C[iO] + C[iR]); Finite(bound, "response bound");
            double nu = bound > 0 ? Math.Min(1, Norm(desired) / bound) : 0;
            var zeta = new double[K]; zeta[iO] = I.chi * nu; zeta[iR] = I.chi * nu + I.lockR;
            return UpdateCChecked(C, I.eta, rho * ell, I.q, validity, I.macc, I.delta0 * (1 - I.phi * rho), zeta, I.OmC);
        }
        static string JsonVector(double[] values) { return "[" + string.Join(",", values.Select(x => x.ToString("R", Inv))) + "]"; }
        static string JsonRows(IEnumerable<double[]> rows) { return "[" + string.Join(",", rows.Select(JsonVector)) + "]"; }
        static string JsonMatrix(double[,] matrix)
        {
            return JsonRows(Enumerable.Range(0, matrix.GetLength(0)).Select(row => Enumerable.Range(0, matrix.GetLength(1)).Select(column => matrix[row, column]).ToArray()));
        }
        public string Export(int id, List<double[]> initial, List<double[]> terminal)
        {
            var text = new StringBuilder("{\"system_id\":" + id + ",\"closure\":\"capacity_subsystem_fixed_forcing\"");
            text.Append(",\"monitoring_form\":\"").Append(I.monitoringForm).Append("\",\"salience_form\":\"").Append(I.salienceForm).Append("\"");
            text.Append(",\"required_routes\":[").Append(string.Join(",", I.reportingRouteRequired.Select(value => value ? "true" : "false"))).Append("]");
            text.Append(",\"available_routes\":[").Append(string.Join(",", I.reportingRouteAvailable.Select(value => value ? "true" : "false"))).Append("]");
            text.Append(",\"external_capacity_change\":").Append(JsonVector(I.OmC ?? new double[K]));
            var scalars = new Dictionary<string, double> {
                {"eta", I.eta}, {"delta0", I.delta0}, {"chi", I.chi}, {"lockR", I.lockR}, {"theta", I.theta}, {"phi", I.phi},
                {"mcov", I.mcov}, {"macc", I.macc}, {"mtime", I.mtime}, {"mplur", I.mplur}, {"gamma", I.gamma}, {"rho", I.rho},
                {"h", I.h}, {"s0", I.s0}, {"kappa", I.kappa}
            };
            foreach (var item in scalars) text.Append(",\"").Append(item.Key).Append("\":").Append(item.Value.ToString("R", Inv));
            text.Append(",\"q\":").Append(JsonVector(I.q)).Append(",\"beta\":").Append(JsonVector(I.beta));
            text.Append(",\"signal\":").Append(JsonVector(signal)).Append(",\"error\":").Append(JsonVector(I.xi));
            text.Append(",\"AI\":").Append(JsonMatrix(I.AI)).Append(",\"AT\":").Append(JsonMatrix(I.AT));
            text.Append(",\"initial_states\":").Append(JsonRows(initial)).Append(",\"terminal_states\":").Append(JsonRows(terminal));
            return text.Append("}").ToString();
        }
    }
    public static double[] StepFeedback(Inp I, double[] C)
    {
        Vector(C, K, "feedback capacity", true, true);
        var result = new FeedbackMap(I).Step(C);
        Vector(result, K, "feedback capacity output", true, true); return result;
    }
    static Inp DrawFeedback(ref int rejected)
    {
        while (true)
        {
            var input = Draw(RB); input.OmC = null; input.delta0 = Math.Max(input.delta0, 0.01);
            if (input.eta + input.delta0 + input.chi + input.lockR <= 1) return input;
            rejected++;
        }
    }
    static void RunDynamics()
    {
        Rng = new Random(20261008); Reset(); double worst = 0;
        for (int t = 0; t < 10000; t++)
        {
            double eta = U(0, 0.9), dl = U(0.001, 1 - eta), zr = U(0, 1 - eta - dl), ellR = U(0, 1), qk = U(0, 1), v = U(0, 1), acc = U(0, 1), c = U(0, 1);
            double a = eta * ellR * qk;
            double bInv = a * (1 - v) * (1 - acc);
            for (int s = 0; s < 20000; s++) c = c + a * v * (1 - c) - bInv * c - (dl + zr) * c;
            double err = Math.Abs(c - a * v / (a * v + bInv + dl + zr)); worst = Math.Max(worst, err);
            Check("T19 Without feedback, capacity converges to the steady state of Property 4", "Property", err <= TolConv);
        }
        // With feedback, convergence and uniqueness are not proven; they are examined numerically from five initial states
        long slow = 0, multi = 0; string multiEx = null; int rejectedFeedback = 0, rejectedTrap = 0;
        for (int t = 0; t < 2000; t++)
        {
            var I = DrawFeedback(ref rejectedFeedback); var map = new FeedbackMap(I);
            var inits = new List<double[]>(); var finals = new List<double[]>();
            for (int st = 0; st < 5; st++) { var C = Vec(K, 0, 1); inits.Add((double[])C.Clone()); for (int s = 0; s < 4000; s++) C = map.Step(C); finals.Add(C); }
            Func<List<double[]>, double> spread = L => { double sp = 0; for (int a1 = 1; a1 < L.Count; a1++) for (int k = 0; k < K; k++) sp = Math.Max(sp, Math.Abs(L[a1][k] - L[0][k])); return sp; };
            bool unique = spread(finals) <= 1e-6;
            if (!unique)
            {
                finals = inits.Select(c0 => { var C = (double[])c0.Clone(); for (int s = 0; s < 200000; s++) C = map.Step(C); return C; }).ToList();
                unique = spread(finals) <= 1e-6;
                if (unique) slow++; else multi++;
            }
            bool fixedPts = finals.All(f => { var g = map.Step(f); return Enumerable.Range(0, K).All(k => Math.Abs(g[k] - f[k]) <= 1e-8); });
            Check("N1 Sampled capacity terminal states satisfy the fixed-point residual tolerance", "Numerical", fixedPts);
            Check("O3 Initial states reach distinct capacity terminal states (not a stability proof)", "Observation", unique);
            if (!unique) Out.AppendLine("FEEDBACK|" + map.Export(t, inits, finals));
            if (!unique && multiEx == null) multiEx = "theta=" + F3(I.theta) + ",indices=" + string.Join("/", finals.Select(f => F3(Dot(I.beta, f))).Distinct());
        }
        for (int t = 0; t < 2000; t++)
        {
            var I = DrawFeedback(ref rejectedTrap); I.theta = 0; var map = new FeedbackMap(I); var C = Vec(K, 0, 1); C[iL] = 0;
            for (int s = 0; s < 4000; s++) C = map.Step(C);
            Check("T20 Capability trap: with theta = 0, zero learning capacity is absorbing (Property 12)", "Property", C.Max() <= 1e-6);
        }
        Dump("DYNAMICS (seed 20261008)", "C");
        Rep("C|worst", worst.ToString("E2", Inv));
        Rep("C|feedback slow/multi", slow + "/" + multi);
        Out.AppendLine("DYNAMIC_DOMAIN|rejected_feedback=" + rejectedFeedback + "|rejected_trap=" + rejectedTrap + "|rule=S29_after_decay_floor");
        Rep("C|bistable example", multiEx);
    }

    // ---------- convexity (Property 13) ----------
    static void RunConvexity()
    {
        Rng = new Random(20261011); Reset();
        for (int t = 0; t < 20000; t++)
        {
            var I = Draw(RB); I.OmC = null; I.xi = new double[p]; var dh = Vec(m, 0.1, 1); double ndh = Norm(dh); for (int i = 0; i < m; i++) dh[i] /= ndh;
            bool convexCase = t % 2 == 0;
            if (convexCase) I.h = U(1.2, 3); else { I.h = 1; I.chi = 0; }
            var sv = Mul(I.M, Mul(I.E, dh)); double y1 = Norm(sv);
            double ystar = I.h > 1 ? I.s0 * Math.Pow((I.h - 1) / (I.h + 1), 1 / I.h) : 3 * I.s0, zmax = ystar / y1;
            Func<double, double> g = z => { var J2 = I.Copy(); J2.d = dh.Select(x => x * z).ToArray(); return Period(J2, false).dC; };
            var R1 = Period(I, false);
            double alpha = I.eta * R1.rhoT * R1.MQ * I.gamma * Dot(I.beta, I.q.Select((qq, k) => qq * (1 - I.C[k])).ToArray());
            var A = new double[r, p]; for (int i = 0; i < r; i++) for (int j = 0; j < p; j++) A[i, j] = I.AI[i, j] + I.C[iT] * I.AT[i, j];
            var eh = sv.Select(x => x / y1).ToArray();
            double b = R1.bound > 0 ? R1.MQ * I.gamma * Norm(Mul(A, eh)) / R1.bound : 0;
            if (convexCase && alpha < I.chi * (I.beta[iO] * I.C[iO] + I.beta[iR] * I.C[iR]) * b) continue;   // sufficient condition of Property 13
            int Ng = 60; double hs = zmax / Ng; bool ok = true;
            for (int i = 1; i < Ng; i++)
            {
                double g0 = g((i - 1) * hs), g1 = g(i * hs), g2 = g((i + 1) * hs), sd = g0 - 2 * g1 + g2, tol = 1e-10 * (Math.Abs(g1) + 1e-6);
                if (convexCase && sd < -tol) ok = false;
                if (!convexCase && sd > tol) ok = false;
            }
            if (convexCase) Check("T21 With h > 1, capacity response is convex in disturbance below the inflection point (Property 13)", "Property", ok);
            else Check("T22 With h = 1 and no option consumption, capacity response is concave in disturbance (Property 13)", "Property", ok);
        }
        Dump("CONVEXITY (seed 20261011)", "V");
    }

    static void RunSelection(int trials)
    {
        Rng = new Random(20261009); Reset();
        for (int t = 0; t < trials; t++)
        {
            int NA = 5; var om = Simplex(5); double tauC = 0.25, Hmax = 1.0;
            var V = new double[NA]; var feas = new bool[NA];
            for (int i = 0; i < NA; i++)
            {
                V[i] = om[0] * U(0, 1) + om[1] * U(0, 1) + om[2] * U(0, 1) + om[3] * U(0, 1) - om[4] * U(0, 1);
                double minCrit = U(0, 1), Hi = U(0, 1.5); bool ir = Rng.NextDouble() < 0.15;
                double exceedance = Hi > Hmax ? 1 : 0;
                var risk = new RiskBounds(exceedance, new[] { exceedance }, exceedance);
                feas[i] = FeasibleAlternative(new[] { minCrit }, new[] { true }, new[] { tauC },
                    new[] { 0.5 }, new[] { 0.0 }, ir, true, risk, 0, new[] { 0.0 }, 0);
            }
            int best = SelectBest(V, feas);
            var ranked = Enumerable.Range(0, NA).Where(i => feas[i]).OrderByDescending(i => V[i]).ToList();
            Check("T23 Selected alternative is feasible (Property 8)", "Property", best < 0 || feas[best]);
            Check("T24 Selection agrees with independent ranking of feasible alternatives", "Implementation", best == (ranked.Count > 0 ? ranked[0] : -1));
            int ub = 0; for (int i = 1; i < NA; i++) if (V[i] > V[ub]) ub = i;
            Check("O1 Unconstrained maximum would be infeasible", "Observation", feas[ub]);
            Check("O2 No feasible alternative exists", "Observation", best >= 0);
        }
        Dump("ALTERNATIVE SELECTION (seed 20261009, " + trials + " decision problems, 5 alternatives each)", "S");
    }

    // ---------- worked example ----------
    public class Design { public string Name; public double mcov, macc, mtime, mplur, gamma, rho; }
    public class ExParams
    {
        public double delta0 = 0.02, eta = 0.6, h = 2, theta = 0.5, phi = 0.5, chi = 0.02, lockR = 0.005, q = 0.5, s0 = 0.4, kappa = 0.25;
        public double sig = 0.8, noise = 0.2, c0 = 0.4, period = 1, saturationCap = 1.5, horizonTolerance = 0.01;
        public double routine = 0, reversibilityFloor = 0.30;
        public string monitoringForm = "product", salienceForm = "hill";
    }
    public static Inp ExampleInp(Design D, ExParams P, double[] C, double[] S, double dScale)
    {
        if (D == null || P == null) throw new ArgumentException("Explicit design and worked-example parameters are required.");
        Vector(C, K, "worked-example capacity", true, true); Vector(S, 2, "worked-example performance", true, true);
        Probability(P.delta0); Probability(P.lockR); Probability(P.c0); Probability(P.horizonTolerance);
        Positive(P.period, "review period"); Nonnegative(P.sig, "annual disturbance"); Nonnegative(dScale, "disturbance multiplier");
        if (double.IsNaN(P.saturationCap) || P.saturationCap <= 0)
            throw new ArgumentOutOfRangeException("saturationCap", "Use a positive cap or positive infinity for uncapped aggregation.");
        var I = new Inp();
        I.S = (double[])S.Clone(); I.w = new[] { 0.2, 0.8 }; I.tau = new[] { 0.45, 0.0 }; I.crit = new[] { true, false };
        I.C = C; I.beta = Enumerable.Repeat(1.0 / K, K).ToArray(); I.q = Enumerable.Repeat(P.q, K).ToArray();
        I.eta = P.eta; I.delta0 = ReviewRate(P.delta0, P.period); I.chi = P.chi; I.lockR = ReviewRate(P.lockR, P.period);
        I.theta = P.theta; I.phi = P.phi; I.h = P.h; I.s0 = P.s0; I.kappa = P.kappa;
        I.d = new[] { Math.Min(P.saturationCap, P.sig * P.period) * dScale }; I.E = new double[,] { { 1 } }; I.M = new double[,] { { 1 }, { 0 } }; I.xi = new[] { 0.0, P.noise };
        I.routine = new[] { P.routine, 0.0 }; I.capacityFloor = new double[K]; I.capacityFloor[iR] = P.reversibilityFloor;
        I.temporaryCritical = new bool[I.S.Length]; I.cumulativeRequired = false;
        I.monitoringForm = P.monitoringForm; I.salienceForm = P.salienceForm;
        I.mcov = D.mcov; I.macc = D.macc; I.mtime = D.mtime; I.mplur = D.mplur; I.gamma = D.gamma; I.rho = D.rho;
        I.AI = new double[,] { { 1, 0 }, { 0, 1 } }; I.AT = new double[,] { { 0.5, 0 }, { 0, 0.5 } };
        I.B = new double[,] { { 0.05 }, { 0 } }; I.G = new double[,] { { 0.2, 0 }, { 0, 0 } }; I.Om = new[] { 0.0, 0.0 }; I.OmC = null;
        I.pi = new[] { 1.0 }; I.Hm = new double[,] { { 0.3 } }; I.Hmax = 10; I.hgMax = new[] { 0.3 }; I.xIrr = 10; I.EQ = 1; I.EC = 1;
        I.designatedGroups = new[] { true }; I.reportingRouteRequired = new[] { true }; I.reportingRouteAvailable = new[] { true };
        ValidateInput(I);
        return I;
    }
    public class LearningThreshold { public double? value; public string status; }
    public static LearningThreshold CriticalLearning(double eta, double Qv, double erosion)
    {
        Probability(eta); Finite(Qv, "Qv"); Nonnegative(erosion, "erosion");
        if (Math.Abs(Qv) > 1) throw new ArgumentOutOfRangeException("Qv");
        if (eta == 0 || Qv <= 0) return new LearningThreshold { value = null, status = "NoPositiveGainAtAnyLearningLevel" };
        double denominator = eta * Qv;
        if (denominator == 0) throw new ArithmeticException("Critical-learning denominator underflowed; higher precision is required.");
        double threshold = erosion / denominator;
        return new LearningThreshold {
            value = double.IsInfinity(threshold) ? (double?)null : threshold,
            status = denominator <= erosion ? "UnattainableWithinUnitLearning" : "FiniteStrictThreshold"
        };
    }
    public class Traj { public double? maxDCd; public double[] Cidx, CidxCf, CidxCf0, SCs, Ssoc, Cfin; public int[] Cls, Stat; public double first_dC, first_dCd, ellR0, ell0, MQ0, v0; public double? crit0; public string criticalLearningStatus; public bool allOk = true, anyMaterial = false; }
    public static Traj RunDesign(Design D, ExParams P, int nPer, Func<int, bool> dist, double[] C0, double[] S0)
    {
        if (nPer <= 0 || dist == null) throw new ArgumentException("A positive review count and disturbance schedule are required.");
        var initial = ExampleInp(D, P, C0, S0, 0);
        var T = new Traj { Cidx = new double[nPer + 1], CidxCf = new double[nPer + 1], CidxCf0 = new double[nPer + 1], SCs = new double[nPer + 1], Ssoc = new double[nPer + 1], Cls = new int[nPer], Stat = new int[nPer] };
        double[] C = (double[])C0.Clone(), Ccf = (double[])C0.Clone(), Ccf0 = (double[])C0.Clone(), S = (double[])S0.Clone(), Scf = (double[])S0.Clone(), Scf0 = (double[])S0.Clone();
        var initialFloors = new double[K]; initialFloors[iR] = P.reversibilityFloor;
        T.allOk = InitialProtectionHolds(initial);
        var beta = Enumerable.Repeat(1.0 / K, K).ToArray();
        T.Cidx[0] = Dot(beta, C); T.CidxCf[0] = T.Cidx[0]; T.CidxCf0[0] = T.Cidx[0]; T.SCs[0] = 0.2 * S[0] + 0.8 * S[1]; T.Ssoc[0] = S[0];
        for (int t = 0; t < nPer; t++)
        {
            bool on = dist(t);
            var I = ExampleInp(D, P, C, S, on ? 1 : 0);
            var A = Period(I, false); var Z = Period(I, true);
            double dCd = Dot(beta, A.Cn) - Dot(beta, Z.Cn), nd = Norm(I.d); if (on && (!T.maxDCd.HasValue || dCd > T.maxDCd.Value)) T.maxDCd = dCd;
            if (t == 0)
            {
                T.first_dC = A.dC; T.first_dCd = dCd; T.ellR0 = A.ellR; T.ell0 = A.ell; T.MQ0 = A.MQ; T.v0 = A.v;
                var qv = new double[K]; var dz = new double[K];
                for (int k = 0; k < K; k++) { qv[k] = I.q[k] * (A.v * (1 - C[k]) - (1 - A.v) * (1 - I.macc) * C[k]); dz[k] = (A.deltaT + A.zeta[k]) * C[k]; }
                var threshold = CriticalLearning(I.eta, Dot(beta, qv), Dot(beta, dz));
                T.crit0 = threshold.value; T.criticalLearningStatus = threshold.status;
            }
            T.Cls[t] = LatentClassOf(A.SB, A.MQ, A.SC, A.SCn, A.SCmin, nd, dCd, A.dC); T.Stat[t] = Status(A.dC);
            if (!ReviewProtectionHolds(A, T.SCs[0], P.horizonTolerance)) T.allOk = false;
            if (nd > EpsD) T.anyMaterial = true;
            C = A.Cn; S = A.Sn;
            var counter = Period(ExampleInp(D, P, Ccf, Scf, 0), false); Ccf = counter.Cn; Scf = counter.Sn;
            var I0 = ExampleInp(D, P, Ccf0, Scf0, 0); I0.xi = new double[I0.xi.Length]; var counter0 = Period(I0, false); Ccf0 = counter0.Cn; Scf0 = counter0.Sn; T.CidxCf0[t + 1] = Dot(beta, Ccf0);
            T.Cidx[t + 1] = Dot(beta, C); T.CidxCf[t + 1] = Dot(beta, Ccf); T.SCs[t + 1] = 0.2 * S[0] + 0.8 * S[1]; T.Ssoc[t + 1] = S[0];
        }
        T.Cfin = C;
        return T;
    }
    public static Design DA = new Design { Name = "A", mcov = 0.5, macc = 0.8, mtime = 0.5, mplur = 0.3, gamma = 0.2, rho = 0.2 };
    public static Design DB = new Design { Name = "B", mcov = 0.9, macc = 0.8, mtime = 0.9, mplur = 0.8, gamma = 0.9, rho = 0.8 };
    static string ClsSeq(int[] c) { return string.Join("", c.Select(x => "XUFGRS"[x])); }

    static void RunExample()
    {
        Out.AppendLine("== WORKED EXAMPLE (latent classes only: X unacceptable, U unverified, F fragile, G capability gain, R robust, S resilient)");
        Out.AppendLine("Empirical capacity verification for every simulation: " + SimulationEvidenceStatus);
        var P = new ExParams(); var C0 = Enumerable.Repeat(P.c0, K).ToArray(); var S0 = new[] { 0.55, 0.60 };
        foreach (var D in new[] { DA, DB })
        {
            var T = RunDesign(D, P, 5, t => true, C0, S0);
            Rep("W|" + D.Name + "|MQ,ell,ellR,v,crit", F4(T.MQ0) + "," + F4(T.ell0) + "," + F4(T.ellR0) + "," + F4(T.v0) + "," + (T.crit0.HasValue ? F4(T.crit0.Value) : T.criticalLearningStatus));
            Rep("W|" + D.Name + "|dC1,dCd1", F4(T.first_dC) + "," + F4(T.first_dCd));
            Rep("W|" + D.Name + "|C5,Ccf5,SC5,social5", F4(T.Cidx[5]) + "," + F4(T.CidxCf[5]) + "," + F4(T.SCs[5]) + "," + F4(T.Ssoc[5]));
            Rep("W|" + D.Name + "|Ccf5noNoise,maxPeriodGain", F4(T.CidxCf0[5]) + "," + (T.maxDCd.HasValue ? F4(T.maxDCd.Value) : "NotAssessed_NoDisturbanceEvent"));
            Rep("W|" + D.Name + "|classes,status", ClsSeq(T.Cls) + "," + string.Join("", T.Stat.Select(s => s > 0 ? "+" : (s < 0 ? "-" : "0"))));
            var L = RunDesign(D, P, 4000, t => true, C0, S0);
            Rep("W|" + D.Name + "|Csteady", F4(L.Cidx[4000]));
        }
        // Jensen gap of the first-period capacity response to a symmetric spread of disturbance (Property 13)
        foreach (var hh in new[] { 1.0, 2.0 })
            foreach (var nz in new[] { 0.0, 0.2 })
            {
                Func<double, double> gJ = dd => { var P2 = new ExParams { h = hh, noise = nz, sig = dd }; return Period(ExampleInp(DB, P2, C0, S0, 1), false).dC; };
                double ystar = hh > 1 ? P.s0 * Math.Pow((hh - 1) / (hh + 1), 1 / hh) : 0;
                double gapIn = 0.5 * (gJ(0.04) + gJ(0.20)) - gJ(0.12), gapOut = 0.5 * (gJ(0.4) + gJ(1.2)) - gJ(0.8);
                Rep("W|jensen|h=" + hh.ToString("F0", Inv) + ",noise=" + nz.ToString("F1", Inv), "ystar=" + F4(ystar) + ",gap(0.12+/-0.08)=" + gapIn.ToString("E2", Inv) + ",gap(0.80+/-0.40)=" + gapOut.ToString("E2", Inv));
            }
    }

    // ---------- sensitivity analysis ----------
    static void RunSensitivity()
    {
        Out.AppendLine("== SENSITIVITY (worked example; five-year horizon; recurring disturbance)");
        var cases = new List<Tuple<string, ExParams>>();
        cases.Add(Tuple.Create("baseline", new ExParams()));
        foreach (var x in new[] { 0.05, 0.08 }) cases.Add(Tuple.Create("delta0=" + x.ToString("F2", Inv), new ExParams { delta0 = x }));
        foreach (var x in new[] { 0.3, 0.9 }) cases.Add(Tuple.Create("eta=" + x.ToString("F1", Inv), new ExParams { eta = x }));
        foreach (var x in new[] { 1.0, 3.0 }) cases.Add(Tuple.Create("h=" + x.ToString("F0", Inv), new ExParams { h = x }));
        foreach (var x in new[] { 1.0, 0.0 }) cases.Add(Tuple.Create("theta=" + x.ToString("F1", Inv), new ExParams { theta = x }));
        foreach (var x in new[] { 0.0, 1.0 }) cases.Add(Tuple.Create("phi=" + x.ToString("F1", Inv), new ExParams { phi = x }));
        foreach (var x in new[] { 0.25, 1.0 }) cases.Add(Tuple.Create("q=" + x.ToString("F2", Inv), new ExParams { q = x }));
        foreach (var x in new[] { 0.05, 0.6 }) cases.Add(Tuple.Create("noise=" + x.ToString("F2", Inv), new ExParams { noise = x }));
        cases.Add(Tuple.Create("chi=0.10", new ExParams { chi = 0.10 }));
        cases.Add(Tuple.Create("period=5y", new ExParams { period = 5 }));
        foreach (var cs in cases)
        {
            var P = cs.Item2; var C0 = Enumerable.Repeat(P.c0, K).ToArray(); var S0 = new[] { 0.55, 0.60 }; int nPer = (int)Math.Round(5 / P.period);
            var TA = RunDesign(DA, P, nPer, t => true, C0, S0); var TB = RunDesign(DB, P, nPer, t => true, C0, S0);
            Rep("Z|" + cs.Item1, "A:C5=" + F3(TA.Cidx[nPer]) + ",cls=" + ClsSeq(TA.Cls) + "|B:C5=" + F3(TB.Cidx[nPer]) + ",gainVsCf=" + F3(TB.Cidx[nPer] - TB.CidxCf[nPer]) + ",cls=" + ClsSeq(TB.Cls) + "|B-A=" + F3(TB.Cidx[nPer] - TA.Cidx[nPer]));
        }
    }

    // ---------- joint (global) sensitivity analysis ----------
    static double[] Ranks(double[] x)
    {
        var idx = Enumerable.Range(0, x.Length).OrderBy(i => x[i]).ToArray(); var r = new double[x.Length];
        for (int i = 0; i < idx.Length; )
        {
            int j = i; while (j + 1 < idx.Length && x[idx[j + 1]] == x[idx[i]]) j++;
            double avg = (i + j) / 2.0 + 1; for (int k = i; k <= j; k++) r[idx[k]] = avg; i = j + 1;
        }
        return r;
    }
    static double Spearman(double[] a, double[] b)
    {
        var ra = Ranks(a); var rb = Ranks(b); double ma = ra.Average(), mb = rb.Average(), sab = 0, saa = 0, sbb = 0;
        for (int i = 0; i < ra.Length; i++) { sab += (ra[i] - ma) * (rb[i] - mb); saa += (ra[i] - ma) * (ra[i] - ma); sbb += (rb[i] - mb) * (rb[i] - mb); }
        return sab / Math.Sqrt(saa * sbb);
    }
    static void RunGlobal(int draws)
    {
        Out.AppendLine("== JOINT SENSITIVITY (" + draws + " joint draws of all illustrative parameters; five-year horizon; seed 20261013)");
        Rng = new Random(20261013);
        string[] names = { "delta0", "eta", "h", "theta", "phi", "q", "noise", "chi", "lock", "sigma0", "kappa", "disturbance", "c0" };
        var X = new double[names.Length][]; for (int k = 0; k < names.Length; k++) X[k] = new double[draws];
        var level = new double[draws]; var gain = new double[draws];
        int bGtA = 0, bGtCf = 0, bAllAf = 0, bAnyAf = 0, bRises = 0, aFav = 0, bGtCf0 = 0; double maxGain = double.NegativeInfinity;
        var S0 = new[] { 0.55, 0.60 };
        Out.AppendLine("GD|draw," + string.Join(",", names) +
            ",A_capacity5,B_capacity5,B_counterfactual5,B_counterfactual_no_error5,A_classes,B_classes,horizon_tolerance,A_protected_horizon,B_protected_horizon,A_retained_gain,B_retained_gain");
        for (int t = 0; t < draws; t++)
        {
            var P = new ExParams { delta0 = U(0.01, 0.10), eta = U(0.2, 0.8), h = U(1, 3), theta = U(0, 1), phi = U(0, 1), q = U(0.1, 1),
                noise = U(0.05, 0.6), chi = U(0, 0.08), lockR = U(0, 0.01), s0 = U(0.2, 0.8), kappa = U(0.1, 0.5), sig = U(0.3, 1.2), c0 = U(0.2, 0.6) };
            double[] v = { P.delta0, P.eta, P.h, P.theta, P.phi, P.q, P.noise, P.chi, P.lockR, P.s0, P.kappa, P.sig, P.c0 };
            for (int k = 0; k < v.Length; k++) X[k][t] = v[k];
            var C0 = Enumerable.Repeat(P.c0, K).ToArray();
            var TA = RunDesign(DA, P, 5, s => true, C0, S0); var TB = RunDesign(DB, P, 5, s => true, C0, S0);
            Out.AppendLine("GD|" + (t + 1).ToString(Inv) + "," +
                string.Join(",", v.Select(value => value.ToString("R", Inv))) + "," +
                string.Join(",", new[] { TA.Cidx[5], TB.Cidx[5], TB.CidxCf[5], TB.CidxCf0[5] }.Select(value => value.ToString("R", Inv))) +
                "," + ClsSeq(TA.Cls) + "," + ClsSeq(TB.Cls) + "," + P.horizonTolerance.ToString("R", Inv) +
                "," + TA.allOk + "," + TB.allOk + "," + RetainedGain(TA) + "," + RetainedGain(TB));
            level[t] = TB.Cidx[5] - TB.Cidx[0]; gain[t] = TB.Cidx[5] - TB.CidxCf[5];
            if (TB.Cidx[5] > TA.Cidx[5]) bGtA++;
            if (gain[t] > 0) bGtCf++;
            if (TB.Cidx[5] > TB.CidxCf0[5]) bGtCf0++;
            if (TB.maxDCd.HasValue && TB.maxDCd.Value > maxGain) maxGain = TB.maxDCd.Value;
            if (TB.Cls.All(c => c == 3)) bAllAf++;
            if (TB.Cls.Any(c => c == 3)) bAnyAf++;
            if (level[t] > 0) bRises++;
            if (TA.Cls.Any(c => c >= 3)) aFav++;
        }
        Func<int, string> pct = c => (100.0 * c / draws).ToString("F1", Inv);
        Rep("G|shares(%)", "BaboveA=" + pct(bGtA) + ",BaboveCounterfactual=" + pct(bGtCf) + ",BcapabilityGainEveryYear=" + pct(bAllAf) + ",BcapabilityGainAnyYear=" + pct(bAnyAf) + ",BlevelRises=" + pct(bRises) + ",AanyFavourable=" + pct(aFav));
        Rep("G|noiseFreeCounterfactual", "BaboveNoiseFreeCf(%)=" + pct(bGtCf0) + ",maxPeriodGainB=" + F4(maxGain));
        foreach (var o in new[] { Tuple.Create("levelChange", level), Tuple.Create("gainVsCounterfactual", gain) })
        {
            var rho = Enumerable.Range(0, names.Length).Select(k => Tuple.Create(names[k], Spearman(X[k], o.Item2))).OrderByDescending(z => Math.Abs(z.Item2)).ToList();
            Rep("G|spearman(" + o.Item1 + ")", string.Join(",", rho.Select(z => z.Item1 + "=" + z.Item2.ToString("F2", Inv))));
        }
    }

    // ---------- horizon analysis with episodic disturbances ----------
    // System-level judgement over a horizon (Eq. S32): safeguards, verification, and maintained performance in every
    // period; at least one material disturbance; retained attributable gain; no erosion of the capacity level.
    public static bool RetainedGain(Traj T)
    {
        if (T == null || T.Cidx == null || T.CidxCf == null || T.Cidx.Length < 2 || T.Cidx.Length != T.CidxCf.Length)
            throw new ArgumentException("Complete matching actual/comparator index histories are required.");
        Vector(T.Cidx, T.Cidx.Length, "actual index history", true, true);
        Vector(T.CidxCf, T.CidxCf.Length, "comparator index history", true, true);
        int end = T.Cidx.Length - 1;
        return T.allOk && T.anyMaterial && T.Cidx[end] - T.CidxCf[end] > EpsC && T.Cidx[end] - T.Cidx[0] >= -EpsC;
    }
    static string Judge(Traj T, int H)
    {
        if (!T.allOk) return "no retained-gain verdict (initial/review protections, monitoring, or fixed-baseline performance failed)";
        if (!T.anyMaterial) return "no retained-gain verdict (no material disturbance)";
        double gain = T.Cidx[H] - T.CidxCf[H], level = T.Cidx[H] - T.Cidx[0];
        if (gain <= EpsC) return "no retained-gain verdict (no retained attributable gain)";
        if (level < -EpsC) return "no retained-gain verdict (capacity index eroded)";
        return "retained capability gain";
    }
    static void RunHorizon()
    {
        Out.AppendLine("== HORIZON (30 annual periods; disturbances in periods 3, 8, 13, 18, 23, 28)");
        var P = new ExParams(); var S0 = new[] { 0.55, 0.60 };
        Func<int, bool> ep = t => t % 5 == 2;
        var C0 = Enumerable.Repeat(P.c0, K).ToArray();
        Func<Traj, string> row = T => F4(T.Cidx[0]) + "," + F4(T.Cidx[30]) + "," + F4(T.CidxCf[30]) + "," + (T.Cidx[30] - T.Cidx[0]).ToString("E2", Inv) + "," + F4(T.Cidx[30] - T.CidxCf[30]) + "," + T.Cls.Count(c => c == 3);
        foreach (var D in new[] { DA, DB })
        {
            var T = RunDesign(D, P, 30, ep, C0, S0);
            Rep("H|" + D.Name + "|C0,C30,Ccf30,levelChange,retainedGain,gainPeriods", row(T));
            Rep("H|" + D.Name + "|classes", ClsSeq(T.Cls));
            Rep("H|" + D.Name + "|judgement", Judge(T, 30));
        }
        // Mature system: start from the periodic state reached under the same episodic disturbances (4,000-period burn-in)
        var Le = RunDesign(DB, P, 4000, ep, C0, S0);
        var Te = RunDesign(DB, P, 30, ep, Le.Cfin, S0);
        Rep("H|B-mature|C0,C30,Ccf30,levelChange,retainedGain,gainPeriods", row(Te));
        Rep("H|B-mature|classes", ClsSeq(Te.Cls));
        Rep("H|B-mature|judgement", Judge(Te, 30));
        var Lb = RunDesign(DB, P, 4000, t => true, C0, S0);
        var M = RunDesign(DB, P, 1, t => true, Lb.Cfin, S0);
        Rep("H|mature|dC,dCd,class", M.first_dC.ToString("E2", Inv) + "," + F4(M.first_dCd) + "," + LatentClassNames[M.Cls[0]]);
    }

    // ---------- deterministic edge cases ----------
    static Inp Base()
    {
        Rng = new Random(20261010);
        var I = new Inp();
        I.S = Enumerable.Repeat(0.5, 2).ToArray(); I.w = new[] { 0.5, 0.5 }; I.tau = new[] { 0.2, 0.0 }; I.crit = new[] { true, false };
        I.C = Vec(K, 0.2, 0.8); I.beta = Simplex(K); I.q = Vec(K, 0.2, 1);
        I.eta = 0.5; I.delta0 = 0.1; I.chi = 0.05; I.lockR = 0.01; I.theta = 0.5; I.phi = 0.5; I.h = 2; I.s0 = 0.5; I.kappa = 1;
        I.d = new[] { 1.0 }; I.E = new double[,] { { 1 } }; I.M = new double[,] { { 1 }, { 0 } }; I.xi = new[] { 0.0, 0.1 };
        I.routine = new double[2]; I.capacityFloor = new double[K]; I.temporaryCritical = new bool[2];
        I.cumulativeRequired = false;
        I.mcov = 1; I.macc = 1; I.mtime = 1; I.mplur = 1; I.gamma = 1; I.rho = 1;
        I.AI = new double[,] { { 1, 0 }, { 0, 1 } }; I.AT = new double[2, 2];
        I.B = new double[,] { { 0 }, { 0 } }; I.G = new double[,] { { 0, 0 }, { 0, 0 } }; I.Om = new[] { 0.0, 0.0 }; I.OmC = null;
        I.pi = new[] { 1.0 }; I.Hm = new double[,] { { 0.1 } }; I.Hmax = 10; I.hgMax = new[] { 0.0 }; I.xIrr = 10; I.EQ = 1; I.EC = 1;
        I.designatedGroups = new[] { false }; I.reportingRouteRequired = new[] { true }; I.reportingRouteAvailable = new[] { true };
        return I;
    }
    static void RunEdgeCases()
    {
        Reset();
        var I = Base(); var b = I.beta;
        { var J1 = I.Copy(); J1.d = new[] { 0.0 }; J1.xi = new[] { 0.0, 0.0 }; var R = Period(J1, false);
          double exp = -Enumerable.Range(0, K).Sum(k => b[k] * (R.deltaT + R.zeta[k]) * I.C[k]);
          Check("E1 Zero evidence: no learning, no response, capacity declines by decay and lock-in", "Edge", R.ell == 0 && Norm(R.u) == 0 && Math.Abs(R.dC - exp) < TolBound); }
        { var J1 = I.Copy(); J1.rho = 0; var R = Period(J1, false); Check("E2 Zero retention: response still possible, no retained learning, no capacity gain", "Edge", R.ell > 0 && Norm(R.u) > 0 && R.ellR == 0 && R.dC <= 0); }
        { var J1 = I.Copy(); J1.gamma = 0; var R = Period(J1, false); Check("E3 Zero connectivity: no learning and no response", "Edge", R.ell == 0 && Norm(R.u) == 0); }
        var ones = Enumerable.Repeat(1.0, K).ToArray(); var zeros = new double[K]; var qz = new double[K];
        Check("E4 Capacities at 1, full valid learning, no decay: remain at 1", "Edge", UpdateC(ones, 1, 1, Enumerable.Repeat(1.0, K).ToArray(), 1, 0, 0, qz, null).All(v => Math.Abs(v - 1) < TolBound));
        Check("E5 Capacities at 0, no learning, no decay: remain at 0", "Edge", UpdateC(zeros, 0.7, 0, I.q, 1, 0, 0, qz, null).All(v => v == 0));
        { var zb = new double[K]; zb[iO] = 0.2; zb[iR] = 0.2; bool ok = true;
          foreach (var C in new[] { ones, zeros, I.C }) foreach (var v in new[] { 0.0, 1.0 }) ok &= UpdateC(C, 0.6, 1, Enumerable.Repeat(1.0, K).ToArray(), v, 0, 0.2, zb, null).All(x => x >= -TolBound && x <= 1 + TolBound);
          Check("E6 Boundary parameters eta + delta + zeta = 1: capacity remains in [0,1]", "Edge", ok); }
        { var J1 = I.Copy(); J1.C = (double[])I.C.Clone(); J1.C[iF] = 0; J1.C[iO] = 0; J1.C[iR] = 0; Check("E7 No flexibility, optionality, or reversibility: no response", "Edge", Norm(Period(J1, false).u) == 0); }
        { var J1 = I.Copy(); J1.d = new[] { 0.0 }; var A = Period(J1, false); var Z = Period(J1, true); double dCd = Dot(b, A.Cn) - Dot(b, Z.Cn);
          Check("E8 No disturbance: no attributable gain, zero index, no capability-gain verdict", "Edge", dCd == 0 && Cgi(dCd, 1, 1, true, 0) == 0 && LatentClassOf(true, 1, 0.5, 0.5, 0.5, 0, dCd, 0) != 3); }
        { var J1 = I.Copy(); J1.d = new[] { 1e6 }; J1.B = new double[,] { { 1 }, { 1 } }; var R = Period(J1, false);
          Check("E9 Extreme disturbance: learning, capacity, and performance remain bounded", "Edge", R.ell <= 1 && R.Cn.All(v => v >= 0 && v <= 1) && R.Sn.All(v => v == 0)); }
        { var J1 = I.Copy(); J1.d = new[] { 0.0 }; J1.xi = new[] { 0.3, 0.4 }; var R = Period(J1, false);
          Check("E10 Pure measurement noise: zero validity and no capacity gain", "Edge", R.v == 0 && R.ell > 0 && R.dC <= 0); }
        Check("E11 Attributable gain with safeguard breach: unacceptable class and zero index", "Edge", LatentClassOf(false, 1, 0.5, 0.6, 0.6, 1, 0.05, 0.05) == 0 && Cgi(0.05, 1, 1, false, 1) == 0);
        Check("E12 Attributable loss with safeguard breach: loss reported in full", "Edge", Math.Abs(Cgi(-0.05, 1, 1, false, 1) - (-0.05 / (1 + Eps0))) < TolBound);
        Check("E13 Gain with monitoring below MQmin: unverified class", "Edge", LatentClassOf(true, MQmin / 2, 0.5, 0.5, 0.3, 1, 0.05, 0.05) == 1);
        Check("E14 Robust: no trough, no attributable gain", "Edge", LatentClassOf(true, 1, 0.5, 0.5, 0.5, 1, 0, 0) == 4);
        Check("E15 Resilient: trough then recovery, no attributable gain", "Edge", LatentClassOf(true, 1, 0.5, 0.5, 0.3, 1, 0, 0) == 5);
        Check("E16 Fragile: no recovery, even with attributable gain", "Edge", LatentClassOf(true, 1, 0.5, 0.3, 0.3, 1, 0.05, 0.05) == 2);
        Check("E17 Capability gain: trough, recovery, and attributable gain", "Edge", LatentClassOf(true, 1, 0.5, 0.5, 0.3, 1, 0.05, 0.05) == 3);
        Check("E18 Steady state: stable capacity level with attributable gain", "Edge", LatentClassOf(true, 1, 0.5, 0.5, 0.5, 1, 0.02, 0) == 3);
        Check("E19 Maintained performance with eroding capacity: robust class, eroding status", "Edge", LatentClassOf(true, 1, 0.5, 0.5, 0.5, 1, 0, -0.03) == 4 && Status(-0.03) == -1);
        Check("E20 Empty feasible set: no alternative selected", "Edge", SelectBest(new[] { 0.9, 0.5, 0.1 }, new[] { false, false, false }) == -1);
        { var J1 = I.Copy(); J1.theta = 0; var C = (double[])I.C.Clone(); C[iL] = 0; for (int s = 0; s < 3000; s++) C = StepFeedback(J1, C);
          Check("E21 Capability trap: theta = 0 and L = 0 lead to zero capacity", "Edge", C.Max() < 1e-9); }
        { var J1 = I.Copy(); J1.OmC = Enumerable.Repeat(-0.03, K).ToArray(); var A = Period(J1, false); var Z = Period(J1, true);
          var A0 = Period(I, false); var Z0 = Period(I, true);
          Check("E22 External capacity change cancels from the attributable gain while truncation is inactive", "Edge", Math.Abs((Dot(b, A.Cn) - Dot(b, Z.Cn)) - (Dot(b, A0.Cn) - Dot(b, Z0.Cn))) < TolBound); }
        Dump("EDGE CASES (deterministic)", "E");
    }

    static void RunInterfaceChecks()
    {
        Reset();
        var established = Enumerable.Repeat(true, 8).ToArray();
        Check("R1 Zero EC cannot satisfy empirical evidence requirements", "Edge",
            CapacityEvidenceStatus(0, false, established) == "InsufficientEvidence");
        Check("R2 Tiny positive EC without the evidence criteria is insufficient", "Edge",
            CapacityEvidenceStatus(0.000001, false, new bool[8]) == "InsufficientEvidence");
        Check("R3 Complete checklist requires independent analysis, not an automatic positive claim", "Edge",
            CapacityEvidenceStatus(1, false, established) == "ReadyForIndependentEffectAndSafeguardAnalysis");
        Check("R4 Simulated states never become empirically verified", "Edge",
            CapacityEvidenceStatus(1, true, established) == SimulationEvidenceStatus);
        var C0 = Enumerable.Repeat(0.4, K).ToArray();
        var S0 = new[] { 0.55, 0.60 };
        var I = ExampleInp(DB, new ExParams(), C0, S0, 1);
        I.EC = 0;
        var A = Period(I, false); var Z = Period(I, true);
        double dCd = Dot(I.beta, A.Cn) - Dot(I.beta, Z.Cn);
        int latent = LatentClassOf(A.SB, A.MQ, A.SC, A.SCn, A.SCmin, Norm(I.d), dCd, A.dC);
        double cgi = Cgi(dCd, I.EC, I.EQ, A.SB, Norm(I.d));
        string evidence = CapacityEvidenceStatus(I.EC, false, established);
        Check("R5 Zero EC latent-gain counterexample remains explicit", "Edge",
            latent == 3 && cgi == 0 && evidence == "InsufficientEvidence");
        Rep("R|EC0|latentClass,CGI,capacityEvidenceStatus",
            LatentClassNames[latent] + "," + F4(cgi) + "," + evidence);
        foreach (bool capped in new[] { true, false })
        {
            var P = new ExParams { period = 5, saturationCap = capped ? 1.5 : double.PositiveInfinity };
            foreach (var D in new[] { DA, DB })
            {
                var J = ExampleInp(D, P, C0, S0, 1);
                var R = Period(J, false); var R0 = Period(J, true);
                double gain = Dot(J.beta, R.Cn) - Dot(J.beta, R0.Cn);
                int cls = LatentClassOf(R.SB, R.MQ, R.SC, R.SCn, R.SCmin, Norm(J.d), gain, R.dC);
                Check("R6 " + (capped ? "Capped" : "Uncapped") + D.Name + " five-year exposure and harm", "Edge",
                    J.d[0] == (capped ? 1.5 : 4.0) && cls == 0 && 0.3 * J.d[0] > 0.3);
                Rep("R|" + (capped ? "cap5" : "uncap5") + "|" + D.Name + "|d,social,harm,C,latentClass",
                    F4(J.d[0]) + "," + F4(R.Sn[0]) + "," + F4(0.3 * J.d[0]) + "," + F4(Dot(J.beta, R.Cn)) + "," + LatentClassNames[cls]);
            }
        }
        Func<double, double, double> compound = (x, y) =>
        {
            var J = ExampleInp(DB, new ExParams { noise = 0 }, C0, S0, 1);
            J.d = new[] { x, y };
            J.E = new double[,] { { 1, 0 }, { 0, 1 } };
            J.M = new double[,] { { 1, 0 }, { 0, 1 } };
            J.AI = new double[,] { { 1, -1 }, { 0, 0 } };
            J.AT = new double[2, 2];
            J.B = new double[,] { { 0.05, 0 }, { 0, 0 } };
            J.Hm = new double[,] { { 0.3, 0 } };
            return Period(J, false).dC;
        };
        double gap = 0.5 * (compound(0.101, 0.099) + compound(0.099, 0.101)) - compound(0.1, 0.1);
        Check("R7 Radial convexity does not imply compound-direction Jensen gain", "Edge",
            gap < 0 && Math.Sqrt(0.101 * 0.101 + 0.099 * 0.099) < 0.4 / Math.Sqrt(3));
        Rep("R|compoundDirection|jensenGap", gap.ToString("E10", Inv));
        Check("R8 One component category step in a six-capacity mean is 1/30", "Edge",
            Math.Abs(((0.6 + 5 * 0.4) / 6 - 0.4) - 1.0 / 30) < TolBound);
        Dump("INTERFACE CHECKS (latent/empirical separation, explicit cap, and domain limits)", "R");
    }

    public class CheckCase
    {
        public string Id, Note;
        public bool Passed;
        public double[] Values;
    }
    static void Case(List<CheckCase> rows, string id, bool passed, double[] values, string note)
    {
        rows.Add(new CheckCase { Id = id, Passed = passed, Values = values, Note = note });
    }
    static bool Near(double a, double b) { return Math.Abs(a - b) <= 1e-12; }
    static bool ArgumentRejected(Action action)
    {
        try { action(); }
        catch (ArgumentException) { return true; }
        return false;
    }
    public static CheckCase[] ScenarioCases()
    {
        var rows = new List<CheckCase>();
        var C = Enumerable.Repeat(0.4, K).ToArray(); var S = new[] { 0.55, 0.60 };
        var beta = Enumerable.Repeat(1.0 / K, K).ToArray();
        var P = new ExParams { routine = 0.2 };
        var routineInput = ExampleInp(DB, P, C, S, 0);
        var routine = Period(routineInput, false);
        Case(rows, "SC01", routine.v > 0 && routine.ellR > 0,
            new[] { routine.v, routine.ellR, Dot(beta, routine.Cn) }, "Valid routine learning remains in the no-disturbance arm.");
        var noExposureSignal = ExampleInp(DB, P, C, S, 1); noExposureSignal.M = new double[2, 1];
        var same = Period(noExposureSignal, false); var same0 = Period(noExposureSignal, true);
        double sameGain = Dot(beta, same.Cn) - Dot(beta, same0.Cn);
        Case(rows, "SC02", Near(sameGain, 0) && same.ellR > 0 &&
            LatentClassOf(same.SB, same.MQ, same.SC, same.SCn, same.SCmin, Norm(noExposureSignal.d), sameGain, same.dC) != 3,
            new[] { sameGain, same.ellR }, "Learning equally in both arms receives no disturbance credit.");
        var actual = Period(ExampleInp(DB, P, C, S, 1), false);
        var zeroRoutine = Period(ExampleInp(DB, new ExParams(), C, S, 0), false);
        Case(rows, "SC03", routine.v > zeroRoutine.v && !Near(Dot(beta, routine.Cn), Dot(beta, zeroRoutine.Cn)),
            new[] { Dot(beta, actual.Cn), Dot(beta, routine.Cn), Dot(beta, actual.Cn) - Dot(beta, routine.Cn),
                    Dot(beta, zeroRoutine.Cn) }, "Routine learning changes the comparator; these are not fitted effects.");
        double[] errors = { 0, 0.1, 0.2, 0.4, 0.6 };
        double paired = 0, independent = 0, pairedSquare = 0, independentSquare = 0;
        foreach (double a in errors)
        {
            var pa = new ExParams { routine = 0.2, noise = a };
            double ca = Dot(beta, Period(ExampleInp(DB, pa, C, S, 1), false).Cn);
            double cz = Dot(beta, Period(ExampleInp(DB, pa, C, S, 0), false).Cn);
            paired += (ca - cz) / errors.Length;
            pairedSquare += (ca - cz) * (ca - cz) / errors.Length;
            foreach (double b in errors)
            {
                double cb = Dot(beta, Period(ExampleInp(DB, new ExParams { routine = 0.2, noise = b }, C, S, 0), false).Cn);
                independent += (ca - cb) / (errors.Length * errors.Length);
                independentSquare += (ca - cb) * (ca - cb) / (errors.Length * errors.Length);
            }
        }
        Case(rows, "SC04", Near(paired, independent),
            new[] { paired, independent, pairedSquare - paired * paired, independentSquare - independent * independent },
            "For this finite error law, coupling changes contrast variance, not its expected value; it does not identify a real counterfactual.");
        var concentrated = DiscreteHarmRisk(new[] { new[] { 1.5, 0.0 } }, new[] { 1.0 },
            new[] { 0.5, 0.5 }, 1, new[] { 1.0, 1.0 }, new[] { true, true });
        bool groupUnsafe = FeasibleAlternative(new[] { 0.6 }, new[] { true }, new[] { 0.4 },
            new[] { 0.5 }, new[] { 0.3 }, false, true, concentrated, 0, new[] { 0.0, 0.0 }, 0);
        Case(rows, "SC05", Near(concentrated.aggregate, 0) && Near(concentrated.group[0], 1) && !groupUnsafe,
            new[] { 0.75, concentrated.aggregate, concentrated.group[0] }, "Aggregate-safe concentrated harm is infeasible.");
        var safeRisk = DiscreteHarmRisk(new[] { new[] { 0.2, 0.2 } }, new[] { 1.0 },
            new[] { 0.5, 0.5 }, 1, new[] { 1.0, 1.0 }, new[] { true, true });
        bool safe = FeasibleAlternative(new[] { 0.6 }, new[] { true }, new[] { 0.4 },
            new[] { 0.5 }, new[] { 0.3 }, false, true, safeRisk, 0, new[] { 0.0, 0.0 }, 0);
        Case(rows, "SC06", SelectBest(new[] { 0.9, 0.8 }, new[] { groupUnsafe, safe }) == 1,
            new[] { (double)SelectBest(new[] { 0.9, 0.8 }, new[] { groupUnsafe, safe }) }, "A higher-value group-unsafe alternative cannot win.");
        Case(rows, "SC07", SelectBest(new[] { 0.9, 0.8 }, new[] { groupUnsafe, groupUnsafe }) == -1,
            new[] { -1.0 }, "No group-safe alternative requires redesign or rejection.");
        var boundary = DiscreteHarmRisk(new[] { new[] { 1.5, 0.0 }, new[] { 0.0, 0.0 } },
            new[] { 0.25, 0.75 }, new[] { 0.5, 0.5 }, 1, new[] { 1.0, 1.0 }, new[] { true, true });
        Case(rows, "SC08", Near(boundary.group[0], 0.25) && RiskAccepts(boundary, 0, new[] { 0.25, 0.0 }, 0.25) &&
            !RiskAccepts(boundary, 0, new[] { 0.125, 0.0 }, 0.25),
            new[] { boundary.aggregate, boundary.group[0], boundary.joint }, "Exact finite-law boundary, not validation of estimated risks.");
        var disjoint = DiscreteHarmRisk(new[] { new[] { 1.5, 0.0 }, new[] { 0.0, 1.5 }, new[] { 0.0, 0.0 } },
            new[] { 0.125, 0.125, 0.75 }, new[] { 0.5, 0.5 }, 1, new[] { 1.0, 1.0 }, new[] { true, true });
        Case(rows, "SC09", Near(disjoint.joint, 0.25) &&
            !RiskAccepts(disjoint, 0, new[] { 0.125, 0.125 }, 0.125) &&
            RiskAccepts(disjoint, 0, new[] { 0.125, 0.125 }, 0.25),
            new[] { disjoint.group[0], disjoint.group[1], disjoint.joint }, "Marginal group budgets do not substitute for the joint budget.");
        var point = new RiskBounds(0.08, new[] { 0.08 }, 0.08);
        var upper = new RiskBounds(0.12, new[] { 0.12 }, 0.12);
        Case(rows, "SC10", RiskAccepts(point, 0.1, new[] { 0.1 }, 0.1) &&
            !RiskAccepts(upper, 0.1, new[] { 0.1 }, 0.1),
            new[] { 0.08, 0.12, 0.1 }, "A point estimate cannot replace a required conservative uncertainty bound.");
        var temporary = ExampleInp(DB, new ExParams(), C, S, 1);
        temporary.B[0, 0] = 0.2; temporary.G[0, 0] = 2;
        var recovered = Period(temporary, false);
        temporary.temporaryCritical[0] = true; var breached = Period(temporary, false);
        Case(rows, "SC11", recovered.Sn[0] >= temporary.tau[0] && recovered.Smin[0] < temporary.tau[0] &&
            recovered.SB && !breached.SB, new[] { recovered.Smin[0], recovered.Sn[0], temporary.tau[0] },
            "Review-date recovery cannot erase a prohibited temporary threshold breach.");
        var cumulative = ExampleInp(DB, new ExParams(), C, S, 1);
        cumulative.cumulativeRequired = true; cumulative.cumulativePrior = new[] { 0.6 };
        cumulative.burdenIncrement = new[] { 0.6 }; cumulative.commonBurden = new[] { 0.0 }; cumulative.cumulativeLimit = new[] { 1.0 };
        var cumulativeResult = Period(cumulative, false);
        Case(rows, "SC12", CumulativeSafe(new[] { 0.0 }, new[] { 0.6 }, new[] { 1.0 }) &&
            !cumulativeResult.SB && Near(cumulativeResult.cumulativeHarm[0], 1.2),
            new[] { 0.6, 0.6, cumulativeResult.cumulativeHarm[0], 1.0 },
            "Period-integrated burdens, not endpoint rates, consume a cumulative budget.");
        var full = new Design { mcov = 1, macc = 1, mtime = 1, mplur = 1, gamma = 1, rho = 1 };
        var erosionP = new ExParams { eta = 0.4, delta0 = 0, chi = 0, lockR = 0.6, theta = 1, phi = 0,
            h = 1, s0 = 1, q = 1, sig = 1, noise = 0, reversibilityFloor = 0.4 };
        var half = Enumerable.Repeat(0.5, K).ToArray();
        var erosion = ExampleInp(full, erosionP, half, S, 1); erosion.B = new double[2, 1]; erosion.Hm = new double[1, 1];
        var erosionResult = Period(erosion, false);
        Case(rows, "SC13", Near(Dot(beta, erosionResult.Cn), 0.55) && Near(erosionResult.Cn[iR], 0.3) && !erosionResult.SB,
            new[] { 0.5, Dot(beta, erosionResult.Cn), 0.5, erosionResult.Cn[iR], 0.4 },
            "Rising mean capacity cannot override a declared critical-capacity floor.");
        erosion.capacityFloor[iR] = 0.25; var permitted = Period(erosion, false); var permitted0 = Period(erosion, true);
        int permittedClass = LatentClassOf(permitted.SB, permitted.MQ, permitted.SC, permitted.SCn, permitted.SCmin, 1,
            Dot(beta, permitted.Cn) - Dot(beta, permitted0.Cn), permitted.dC);
        Case(rows, "SC14", permitted.SB && permittedClass == 3,
            new[] { permitted.Cn[iR], 0.25, (double)permittedClass },
            "Above the predeclared floor, component trade-offs remain allowed; not every capacity is required to increase.");
        Case(rows, "SC15", FloorsHold(new[] { 0.4 }, new[] { 0.4 }) && !FloorsHold(new[] { 0.3999 }, new[] { 0.4 }),
            new[] { 0.4, 0.3999 }, "Capacity-floor equality is accepted and a smaller value rejected.");
        var original = CanonicalEvidence(new[] { 20.0 }, new[] { 100.0 }, new[] { 0 }, new[] { 1.0 }, new[] { "source-a:v1" });
        var rescaled = CanonicalEvidence(new[] { 2000.0 }, new[] { 10000.0 }, new[] { 0 }, new[] { 1.0 }, new[] { "source-a:v1" });
        var repeated = CanonicalEvidence(new[] { 20.0, 20.0 }, new[] { 100.0, 100.0 }, new[] { 0, 0 }, new[] { 1.0, 1.0 }, new[] { "source-a:v1", "source-a:v1" });
        Case(rows, "SC16", Near(original[0], rescaled[0]) && Near(original[0], repeated[0]),
            new[] { original[0], rescaled[0], repeated[0] }, "Unit changes and redundant records preserve canonical evidence.");
        Case(rows, "SC17", Salience(Norm(new[] { 0.2, 0.2 }), 0.4, 2) > Salience(Norm(original), 0.4, 2),
            new[] { Salience(Norm(original), 0.4, 2), Salience(Norm(new[] { 0.2, 0.2 }), 0.4, 2) },
            "Unfused duplication deliberately fails invariance; independent canonical dimensions are not interchangeable.");
        Case(rows, "SC18", ArgumentRejected(() => CanonicalEvidence(new[] { 1.0 }, new[] { 0.0 }, new[] { 0 }, new[] { 1.0 }, new[] { "source-a:v1" })) &&
            ArgumentRejected(() => CanonicalEvidence(new[] { 1.0 }, new[] { 1.0 }, new[] { 1 }, new[] { 1.0 }, new[] { "source-a:v1" })),
            new double[0], "Invalid reference units and missing canonical-channel declarations are explicit errors.");
        double[] monitored = { MonitoringQuality(0, 0.9, 0.8, "product"), MonitoringQuality(0, 0.9, 0.8, "geometric"),
            MonitoringQuality(0, 0.9, 0.8, "minimum") };
        Case(rows, "SC19", monitored.All(x => x == 0), monitored, "Every tested monitoring form preserves a missing-channel veto.");
        foreach (string form in new[] { "product", "geometric", "minimum" })
        {
            var formInput = ExampleInp(DB, new ExParams { routine = 0.2, monitoringForm = form }, C, S, 1);
            var f = Period(formInput, false); var f0 = Period(formInput, true);
            Case(rows, "SC20_" + form, f.MQ >= 0 && f.MQ <= 1 && f.Cn.All(x => x >= 0 && x <= 1),
                new[] { f.MQ, Dot(beta, f.Cn), Dot(beta, f0.Cn), Dot(beta, f.Cn) - Dot(beta, f0.Cn) },
                "Structural comparison: gains and magnitudes depend on the aggregation form.");
        }
        foreach (string form in new[] { "hill", "linear_capped" })
        {
            var formInput = ExampleInp(DB, new ExParams { routine = 0.2, salienceForm = form }, C, S, 1);
            var f = Period(formInput, false); var f0 = Period(formInput, true);
            Case(rows, "SC21_" + form, f.sigma >= 0 && f.sigma <= 1,
                new[] { f.sigma, Dot(beta, f.Cn), Dot(beta, f0.Cn), Dot(beta, f.Cn) - Dot(beta, f0.Cn) },
                "The convexity theorem is restricted to its Hill specification and stated conditions.");
        }
        var funded = ExampleInp(DB, new ExParams { routine = 0.2 }, C, S, 1);
        var fundingResult = Period(funded, false);
        funded.gamma *= 0.25; funded.rho *= 0.25; var constrained = Period(funded, false);
        Case(rows, "SC22", Dot(beta, constrained.Cn) < Dot(beta, fundingResult.Cn),
            new[] { Dot(beta, fundingResult.Cn), Dot(beta, constrained.Cn) },
            "A coupled resource shortfall reduces connectivity and retention in this specified scenario.");
        var transformative = ExampleInp(DB, new ExParams(), C, S, 1);
        var proportional = Period(transformative, false);
        transformative.AT = new double[,] { { 0, 1 }, { -1, 0 } }; var distinct = Period(transformative, false);
        double cross = proportional.u[0] * distinct.u[1] - proportional.u[1] * distinct.u[0];
        Case(rows, "SC23", Math.Abs(cross) > 1e-12 && Norm(distinct.u) <= distinct.bound + TolBound,
            new[] { proportional.u[0], proportional.u[1], distinct.u[0], distinct.u[1], distinct.bound },
            "Exogenously available non-proportional response options are not an endogenous theory of institutional transformation.");
        var relevance = Enumerable.Repeat(1.0, K).ToArray();
        double trueChange = Dot(beta, UpdateC(C, 0.6, 0.4, relevance, 0, 0, 0, new double[K], null)) - 0.4;
        double proxyChange = Dot(beta, UpdateC(C, 0.6, 0.4, relevance, 1, 0, 0, new double[K], null)) - 0.4;
        Case(rows, "SC24", trueChange < 0 && proxyChange > 0, new[] { trueChange, proxyChange },
            "Misclassified validity can reverse a predicted change; observed capability contrasts require independent measurement.");
        bool invalid = ArgumentRejected(() => RiskAccepts(new RiskBounds(double.NaN, new[] { 0.0 }, 0), 0, new[] { 0.0 }, 0)) &&
            ArgumentRejected(() => CumulativeSafe(new[] { 0.0 }, new[] { -0.1 }, new[] { 1.0 })) &&
            ArgumentRejected(() => FloorsHold(new[] { 0.5 }, new double[0]));
        Case(rows, "SC25", invalid, new double[0], "Invalid risks, negative burdens and incomplete floor policies cannot silently pass.");
        var random = new Random(20261023); int concentratedCount = 0, emptyCount = 0; bool rankingOk = true;
        for (int t = 0; t < 1000; t++)
        {
            var values = new double[2]; var flags = new bool[2]; int oracle = -1;
            for (int i = 0; i < 2; i++)
            {
                double a = 2 * random.NextDouble(), b = 2 * random.NextDouble(); values[i] = random.NextDouble();
                var risk = DiscreteHarmRisk(new[] { new[] { a, b } }, new[] { 1.0 },
                    new[] { 0.5, 0.5 }, 1, new[] { 1.0, 1.0 }, new[] { true, true });
                flags[i] = FeasibleAlternative(new[] { 0.6 }, new[] { true }, new[] { 0.4 },
                    new[] { 0.5 }, new[] { 0.3 }, false, true, risk, 0, new[] { 0.0, 0.0 }, 0);
                bool direct = a <= 1 && b <= 1 && 0.5 * (a + b) <= 1;
                rankingOk &= flags[i] == direct;
                if (direct && (oracle < 0 || values[i] > values[oracle])) oracle = i;
                if (0.5 * (a + b) <= 1 && (a > 1 || b > 1)) concentratedCount++;
            }
            int selected = SelectBest(values, flags); rankingOk &= selected == oracle;
            if (selected < 0) emptyCount++;
        }
        Case(rows, "SC26", rankingOk && concentratedCount > 0 && emptyCount > 0,
            new[] { 1000.0, (double)concentratedCount, (double)emptyCount },
            "Two-group selection agrees with direct inequalities and independent ranking on the specified seeded test set.");
        bool pathRejected = !FeasibleAlternative(new[] { 0.6 }, new[] { true }, new[] { 0.4 },
            new[] { 0.5 }, new[] { 0.3 }, false, CumulativeSafe(new[] { 0.6 }, new[] { 0.6 }, new[] { 1.0 }),
            safeRisk, 0, new[] { 0.0, 0.0 }, 0);
        Case(rows, "SC27", pathRejected, new[] { 1.2, 1.0 },
            "The alternative selector also vetoes a prohibited cumulative or temporary path even when endpoint and harm-risk tests pass.");
        cumulative.commonBurden = new[] { 0.2 }; cumulative.burdenIncrement = new[] { 0.1 };
        var burdenActual = Period(cumulative, false); var burdenComparator = Period(cumulative, true);
        Case(rows, "SC28", Near(burdenActual.cumulativeHarm[0], 0.9) && Near(burdenComparator.cumulativeHarm[0], 0.8),
            new[] { burdenActual.cumulativeHarm[0], burdenComparator.cumulativeHarm[0] },
            "Common integrated burdens remain in both arms; only the disturbance-linked increment is removed.");
        return rows.ToArray();
    }

    public class AccuracyAssessment { public double? value; public string status; }
    public static AccuracyAssessment InvalidLessonSensitivity(int rejectedInvalid, int retainedInvalid)
    {
        if (rejectedInvalid < 0 || retainedInvalid < 0) throw new ArgumentOutOfRangeException("independently adjudicated lesson counts");
        long total = (long)rejectedInvalid + retainedInvalid;
        return total == 0
            ? new AccuracyAssessment { value = null, status = "Unassessed_NoIndependentlyInvalidLessons" }
            : new AccuracyAssessment { value = rejectedInvalid / (double)total, status = "Assessed_InvalidLessonRejectionSensitivity" };
    }
    public static CheckCase[] SafeguardCases()
    {
        var rows = new List<CheckCase>();
        var c = Enumerable.Repeat(0.4, K).ToArray(); var s = new[] { 0.55, 0.60 };
        var mixed = CanonicalEvidence(new[] { 0.2, 0.6 }, new[] { 1.0, 1.0 }, new[] { 0, 0 }, new[] { 1.0, 1.0 }, new[] { "a:v1", "b:v1" });
        var copied = CanonicalEvidence(new[] { 0.2, 0.6, 0.2 }, new[] { 1.0, 1.0, 1.0 }, new[] { 0, 0, 0 }, new[] { 1.0, 1.0, 1.0 }, new[] { "a:v1", "b:v1", "a:v1" });
        Case(rows, "AC01_mixed_source_copy", Near(mixed[0], 0.4) && Near(copied[0], 0.4) && Near(Salience(Norm(mixed), 0.4, 2), Salience(Norm(copied), 0.4, 2)), new[] { mixed[0], copied[0] }, "Copies do not increase an original source's fixed contribution.");
        Case(rows, "AC02_source_conflict", ArgumentRejected(() => CanonicalEvidence(new[] { 0.2, 0.6 }, new[] { 1.0, 1.0 }, new[] { 0, 0 }, new[] { 1.0, 1.0 }, new[] { "a:v1", "a:v1" })), new double[0], "Conflicting source versions require reconciliation.");
        Case(rows, "AC03_large_norm", Near(Norm(new[] { 1e200 }) / 1e200, 1), new[] { Norm(new[] { 1e200 }) }, "Scaled norm avoids squaring overflow.");
        Case(rows, "AC04_tiny_hill", Salience(1e-200, 1e-200, 2) == 0.5, new[] { Salience(1e-200, 1e-200, 2) }, "Equal tiny strengths and scales give one half.");
        Case(rows, "AC05_large_hill", Salience(1e200, 1, 2) == 1, new[] { Salience(1e200, 1, 2) }, "Log-ratio Hill remains finite at large strength.");
        Case(rows, "AC06_large_validity", Validity(new[] { 1e308 }, new[] { 1e308 }) == 0.5, new[] { Validity(new[] { 1e308 }, new[] { 1e308 }) }, "Max-scaled validity does not overflow the denominator.");
        Case(rows, "AC07_risk_one_rejected", ArgumentRejected(() => RiskAccepts(new RiskBounds(1, new[] { 1.0 }, 1), 1, new[] { 1.0 }, 1)), new double[0], "Upper bounds may equal one, risk tolerances may not.");
        Case(rows, "AC08_risk_bound_one_valid", !RiskAccepts(new RiskBounds(1, new[] { 1.0 }, 1), 0.99, new[] { 0.99 }, 0.99), new[] { 1.0 }, "A valid bound of one fails a smaller valid tolerance.");
        Case(rows, "AC09_eta_zero", CriticalLearning(0, 0.5, 0).value == null && CriticalLearning(0, 0.5, 0).status == "NoPositiveGainAtAnyLearningLevel", new double[0], "Do not divide by eta=0.");
        Case(rows, "AC10_nonpositive_Qv", CriticalLearning(0.6, 0, 0.1).value == null && CriticalLearning(0.6, -0.2, 0.1).value == null, new double[0], "No positive gain when eta Qv is nonpositive and external change is zero.");
        Case(rows, "AC11_unattainable_threshold", CriticalLearning(0.1, 0.2, 0.1).status == "UnattainableWithinUnitLearning", new[] { CriticalLearning(0.1, 0.2, 0.1).value.Value }, "A threshold above one is not attainable retained learning.");
        var accuracy = InvalidLessonSensitivity(20, 0);
        Case(rows, "AC12_accuracy_denominator", accuracy.value == 1, new[] { accuracy.value.Value, 20.0 / 100 }, "Twenty rejected among twenty independently invalid lessons gives sensitivity one, not 0.2.");
        Case(rows, "AC13_accuracy_unassessed", InvalidLessonSensitivity(0, 0).value == null, new double[0], "No invalid adjudicated lessons leaves sensitivity unassessed.");
        var initialBelow = RunDesign(DB, new ExParams(), 5, t => true, c, new[] { 0.449, 0.60 });
        Case(rows, "AC14_initial_performance", !initialBelow.allOk && initialBelow.Cls.All(x => x == 3), new[] { initialBelow.Ssoc[0], initialBelow.Ssoc[1] }, "Later recovery does not retroactively protect the initial horizon.");
        var ratchet = ExampleInp(DB, new ExParams(), c, new[] { 0.85, 0.60 }, 1);
        ratchet.B = new double[2, 1]; ratchet.G = new double[2, 2]; ratchet.Hm = new double[1, 1]; ratchet.Om = new[] { -0.045, 0.0 };
        double initialSC = Dot(ratchet.w, ratchet.S); bool local = true, horizon = true;
        for (int t = 0; t < 5; t++)
        {
            var a = Period(ratchet, false); var z = Period(ratchet, true);
            local &= LatentClassOf(a.SB, a.MQ, a.SC, a.SCn, a.SCmin, Norm(ratchet.d), Dot(ratchet.beta, a.Cn) - Dot(ratchet.beta, z.Cn), a.dC) == 3;
            horizon &= ReviewProtectionHolds(a, initialSC, 0.01); ratchet.C = a.Cn; ratchet.S = a.Sn;
        }
        Case(rows, "AC15_fixed_baseline_ratchet", local && !horizon && Near(Dot(ratchet.w, ratchet.S), 0.605), new[] { initialSC, Dot(ratchet.w, ratchet.S) }, "Five locally tolerated losses cannot pass fixed-baseline horizon protection.");
        var equality = new Res { SB = true, MQ = 0.8, SC = 0.65, SCn = 0.64, SCmin = 0.64, dC = 0 };
        Case(rows, "AC16_horizon_equality", ReviewProtectionHolds(equality, 0.65, 0.01), new[] { 0.65, 0.64, 0.01 }, "Equality at the predeclared fixed-baseline tolerance is accepted.");
        var missingRoute = ExampleInp(DB, new ExParams(), c, s, 1);
        missingRoute.pi = new[] { 1.0, 1.0 }; missingRoute.Hm = new double[,] { { 0.1 }, { 0.1 } };
        missingRoute.hgMax = new[] { 1.0, 1.0 }; missingRoute.designatedGroups = new[] { true, true };
        missingRoute.reportingRouteRequired = new[] { true, true }; missingRoute.reportingRouteAvailable = new[] { true, false };
        foreach (string form in new[] { "product", "geometric", "minimum" })
        {
            missingRoute.monitoringForm = form; var a = Period(missingRoute, false);
            Case(rows, "AC17_required_route_" + form, a.SB && a.MQ == 0 && a.ell == 0 && !a.requiredRoutesPass, new[] { a.MQ, a.ell }, "A missing required group route vetoes aggregate monitoring quality.");
        }
        var overshoot = ExampleInp(DB, new ExParams(), c, s, 1); overshoot.Om = new[] { 1.0, 0.0 }; overshoot.OmC = Enumerable.Repeat(0.9, K).ToArray();
        var raw = Period(overshoot, false);
        Case(rows, "AC18_raw_projection", raw.rawPerformance[0] > 1 && raw.Sn[0] == 1 && raw.rawCapacity[0] > 1 && raw.Cn[0] == 1, new[] { raw.rawPerformance[0], raw.Sn[0], raw.rawCapacity[0], raw.Cn[0] }, "Export raw overshoots; projection is not physical path-safety evidence.");
        Case(rows, "AC19_generic_split_utility", UpdateC(new double[7], 0.6, 0.5, Enumerable.Repeat(0.5, 7).ToArray(), 1, 1, 0.02, new double[7], null).Length == 7, new[] { 7.0 }, "Generic capacity utility still supports the split-component property.");
        var mutations = new Dictionary<string, Action<Inp>> {
            {"S_dimension", i => i.S = new double[0]}, {"S_domain", i => i.S = new[] { -0.1, 0.6 }},
            {"weights_sum", i => i.w = new[] { 0.2, 0.2 }}, {"weights_negative", i => i.w = new[] { -0.2, 1.2 }},
            {"threshold_nonfinite", i => i.tau = new[] { 0.45, double.NaN }}, {"critical_dimension", i => i.crit = new bool[1]},
            {"capacity_dimension", i => i.C = new double[5]}, {"capacity_nonfinite", i => i.C = Enumerable.Repeat(double.NaN, K).ToArray()},
            {"beta_sum", i => i.beta = new double[K]}, {"q_domain", i => i.q = Enumerable.Repeat(1.1, K).ToArray()},
            {"eta_negative", i => i.eta = -0.1}, {"S29", i => i.eta = 1}, {"theta_domain", i => i.theta = 1.1},
            {"phi_domain", i => i.phi = -0.1}, {"rho_nonfinite", i => i.rho = double.PositiveInfinity},
            {"monitoring_domain", i => i.mcov = 1.1}, {"accuracy_domain", i => i.macc = -0.1},
            {"Hill_domain", i => i.h = 0.5}, {"scale_zero", i => i.s0 = 0}, {"response_scale_zero", i => i.kappa = 0},
            {"disturbance_negative", i => i.d = new[] { -0.1 }}, {"exposure_negative", i => i.E = new double[,] { { -1 } }},
            {"monitoring_map_dimension", i => i.M = new double[1, 1]}, {"error_dimension", i => i.xi = new double[1]},
            {"routine_missing", i => i.routine = null}, {"AI_dimension", i => i.AI = new double[1, 1]},
            {"AT_nonfinite", i => i.AT = new double[,] { { double.NaN, 0 }, { 0, 0 } }},
            {"loss_negative", i => i.B = new double[,] { { -1 }, { 0 } }}, {"gain_dimension", i => i.G = new double[1, 1]},
            {"external_performance", i => i.Om = new[] { double.NaN, 0.0 }}, {"external_capacity_dimension", i => i.OmC = new double[1]},
            {"equity_negative", i => i.pi = new[] { -1.0 }}, {"harm_map_negative", i => i.Hm = new double[,] { { -1 } }},
            {"harm_limit_nonfinite", i => i.Hmax = double.NaN}, {"group_limit_nonfinite", i => i.hgMax = new[] { double.NaN }},
            {"group_mask_missing", i => i.designatedGroups = null}, {"route_policy_missing", i => i.reportingRouteRequired = null},
            {"temporary_policy_missing", i => i.temporaryCritical = null}, {"cumulative_required_missing", i => i.cumulativeRequired = true},
            {"irreversibility_nonfinite", i => i.xIrr = double.PositiveInfinity}, {"EC_nonfinite", i => i.EC = double.NaN}
        };
        foreach (var mutation in mutations)
        {
            var input = ExampleInp(DB, new ExParams(), (double[])c.Clone(), s, 1); mutation.Value(input);
            Case(rows, "AC20_" + mutation.Key, ArgumentRejected(() => Period(input, false)), new double[0], "Invalid public inputs must raise a domain/dimension error before labels.");
        }
        foreach (double invalid in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            for (int field = 0; field < 7; field++)
            {
                var a = new[] { 0.8, 0.5, 0.5, 0.5, 0.8, 0.1, 0.0 }; a[field] = invalid;
                Case(rows, "AC21_classifier_" + field + "_" + (double.IsNaN(invalid) ? "NaN" : invalid > 0 ? "PosInf" : "NegInf"),
                    ArgumentRejected(() => LatentClassOf(true, a[0], a[1], a[2], a[3], a[4], a[5], a[6])), new double[0], "No nonfinite scalar may enter a favourable classification.");
            }
        Case(rows, "AC22_dot_dimension", ArgumentRejected(() => Dot(new[] { 1.0 }, new[] { 1.0, 100 })), new double[0], "A trailing dot-product entry cannot be ignored.");
        Case(rows, "AC23_update_domain", ArgumentRejected(() => UpdateC(c, 4, 1, Enumerable.Repeat(1.0, K).ToArray(), 1, 1, 0, new double[K], null)), new double[0], "The boundedness theorem is not applied outside its parameter domain.");
        Case(rows, "AC24_converted_S29", ArgumentRejected(() => ExampleInp(DB, new ExParams { eta = 0.95, period = 5 }, c, s, 1)), new double[0], "Annual admissibility does not replace review-interval admissibility.");
        double largest = 0; var random = new Random(20261023);
        for (int t = 0; t < 30; t++)
        {
            var state = Enumerable.Range(0, K).Select(k => random.NextDouble()).ToArray();
            var input = ExampleInp(DB, new ExParams(), state, s, 1); input.AT = new double[,] { { 0, 1 }, { -1, 0 } };
            var full = Period(input, false).Cn; var reduced = StepFeedback(input, state);
            for (int k = 0; k < K; k++) largest = Math.Max(largest, Math.Abs(full[k] - reduced[k]));
        }
        Case(rows, "AC25_feedback_map_parity", largest <= 1e-12, new[] { largest }, "The cached fixed-forcing capacity map agrees with the full period implementation.");
        return rows.ToArray();
    }

    public static CheckCase[] BoundaryCases()
    {
        var rows = new List<CheckCase>();
        var zero = new double[K]; var ones = Enumerable.Repeat(1.0, K).ToArray();
        var stock = Enumerable.Repeat(0.3, K).ToArray();
        var next = UpdateC(stock, 0.9, 1, ones, 0, 0, 0.1, zero, null);
        Case(rows, "FR01_capacity_boundary_utility", next.All(x => x == 0), next,
            "The admitted zero-survival boundary returns exact zero without stock projection.");
        var design = new Design { mcov = 1, macc = 0, mtime = 1, mplur = 1, gamma = 1, rho = 1 };
        var parameters = new ExParams { eta = 0.9, delta0 = 0.1, lockR = 0, chi = 0,
            phi = 0, theta = 1, q = 1, sig = 0, noise = 1, reversibilityFloor = 0,
            salienceForm = "linear_capped" };
        var boundary = Period(ExampleInp(design, parameters, stock, new[] { 0.55, 0.60 }, 0), false);
        Case(rows, "FR02_capacity_boundary_period", boundary.Cn.All(x => x == 0)
            && boundary.rawCapacity.All(x => x == 0), boundary.Cn, "Full-period/raw parity at the same boundary.");

        double[] values = { 1e308, 1, -1e308 }; string[] identities = { "a:v1", "b:v1", "c:v1" };
        bool invariant = true; int fusions = 0;
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) for (int c = 0; c < 3; c++)
        {
            if (a == b || a == c || b == c) continue;
            var order = new[] { a, b, c };
            double original = CanonicalEvidence(order.Select(i => values[i]).ToArray(), new[] { 1.0, 1.0, 1.0 },
                new[] { 0, 0, 0 }, new[] { 1.0, 1.0, 1.0 }, order.Select(i => identities[i]).ToArray())[0];
            invariant &= original == 1.0 / 3; fusions++;
            for (int copied = 0; copied < 3; copied++) for (int position = 0; position < 4; position++)
            {
                var records = order.ToList(); records.Insert(position, copied);
                double fused = CanonicalEvidence(records.Select(i => values[i]).ToArray(),
                    new[] { 1.0, 1.0, 1.0, 1.0 }, new[] { 0, 0, 0, 0 },
                    new[] { 1.0, 1.0, 1.0, 1.0 }, records.Select(i => identities[i]).ToArray())[0];
                invariant &= fused == 1.0 / 3; fusions++;
            }
        }
        Case(rows, "FR03_signed_fusion_invariance", invariant && fusions == 78,
            new[] { 1.0 / 3, (double)fusions }, "All six permutations and 72 duplicate insertions preserve the finite residual.");
        double longRate = ReviewRate(1e-18, 1e18), shortRate = ReviewRate(0.2, 1e-18);
        Case(rows, "FR04_tiny_annual_rate", Math.Abs(longRate - 0.6321205588285577) <= 2e-16,
            new[] { longRate }, "Tiny annual rates remain effective over long review intervals.");
        Case(rows, "FR05_tiny_review_interval", Math.Abs(shortRate - 2.231435513142098e-19) <= 1e-33,
            new[] { shortRate }, "Short-interval conversion does not collapse to zero.");
        Case(rows, "FR06_zero_rate", ReviewRate(0, 1e308) == 0, new[] { ReviewRate(0, 1e308) }, "Exact zero endpoint.");
        Case(rows, "FR07_unit_rate", ReviewRate(1, 1e-308) == 1, new[] { ReviewRate(1, 1e-308) }, "Exact unit endpoint.");
        Case(rows, "FR08_one_year_identity", ReviewRate(0.1, 1) == 0.1, new[] { ReviewRate(0.1, 1) }, "One-year conversion preserves the input.");
        Case(rows, "FR09_converted_rate_guard", ArgumentRejected(() => ExampleInp(DB,
            new ExParams { delta0 = 1e-18, period = 1e18, eta = 0.7, lockR = 0, chi = 0,
                phi = 0, sig = 0, noise = 0 }, stock, new[] { 0.55, 0.60 }, 0)),
            new[] { 0.7 + longRate }, "Admissibility uses the accurately converted review rate.");

        double[] probabilities = { 0.4285714285714286, 0.28571428571428575, 0.28571428571428575 };
        double belowOne = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(1.0) - 1);
        var certain = DiscreteHarmRisk(new[] { new[] { 1.0 }, new[] { 1.0 }, new[] { 1.0 } },
            probabilities, new[] { 1.0 }, 0, new[] { 0.0 }, new[] { true });
        Case(rows, "FR10_certain_harm", certain.aggregate == 1 && certain.group[0] == 1 && certain.joint == 1
            && !RiskAccepts(certain, belowOne, new[] { belowOne }, belowOne),
            new[] { certain.aggregate, certain.group[0], certain.joint, belowOne }, "Certain harm rejects every tolerance below one.");
        var withZero = DiscreteHarmRisk(new[] { new[] { 1.0 }, new[] { 1.0 }, new[] { 1.0 }, new[] { 0.0 } },
            probabilities.Concat(new[] { 0.0 }).ToArray(), new[] { 1.0 }, 0, new[] { 0.0 }, new[] { true });
        Case(rows, "FR11_zero_probability_safe_scenario", withZero.aggregate == 1 && withZero.group[0] == 1
            && withZero.joint == 1, new[] { withZero.aggregate, withZero.group[0], withZero.joint },
            "A zero-probability safe scenario cannot destroy the certain-harm endpoint.");
        var safe = DiscreteHarmRisk(new[] { new[] { 0.0 }, new[] { 0.0 }, new[] { 0.0 } },
            probabilities, new[] { 1.0 }, 0, new[] { 0.0 }, new[] { true });
        Case(rows, "FR12_zero_harm", safe.aggregate == 0 && safe.group[0] == 0 && safe.joint == 0,
            new[] { safe.aggregate, safe.group[0], safe.joint }, "Exact no-event endpoint.");
        double adjacent = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(1e200) + 1);
        double hill = Salience(adjacent, 1e200, 1e16);
        Case(rows, "FR13_near_threshold_hill", Math.Abs(hill - 0.8454879170168209) <= 2e-15,
            new[] { hill }, "Adjacent strengths agree with an independent high-precision reference.");
        double near1024 = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(1024.0) + 1);
        double nearOne = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(1.0) + 1);
        double originalHill = Salience(near1024, 1024, 1e16), rescaledHill = Salience(nearOne, 1, 1e16);
        Case(rows, "FR14_hill_unit_invariance", originalHill == rescaledHill,
            new[] { originalHill, rescaledHill }, "Exact power-of-two unit rescaling preserves the log ratio.");
        Case(rows, "FR15_above_floor_tradeoff", CriticalHold(new[] { 0.8, 0.4 }, new[] { true, false }, new[] { 0.4, 0.0 })
            && CriticalHold(new[] { 0.7, 0.5 }, new[] { true, false }, new[] { 0.4, 0.0 })
            && Near(0.5 * (0.8 + 0.4), 0.5 * (0.7 + 0.5)), new[] { 0.6, 0.6 },
            "A binding critical floor does not prohibit an above-floor tradeoff.");
        Case(rows, "FR16_critical_breach_veto", !CriticalHold(new[] { 0.39, 1.0 }, new[] { true, false }, new[] { 0.4, 0.0 }),
            new[] { 0.695 }, "A higher aggregate cannot compensate for a critical-floor breach.");
        var low = Enumerable.Repeat(0.1, K).ToArray();
        double clean = UpdateC(low, 0.5, Salience(0.1, 1, 2), ones, 1, 0.8, 0, zero, null)[0] - 0.1;
        double noisy = UpdateC(low, 0.5, Salience(Norm(new[] { 0.1, 0.5 }), 1, 2), ones,
            Validity(new[] { 0.1, 0.0 }, new[] { 0.0, 0.5 }), 0.8, 0, zero, null)[0] - 0.1;
        double pure = UpdateC(low, 0.5, Salience(0.5, 1, 2), ones, 0, 0.8, 0, zero, null)[0] - 0.1;
        Case(rows, "FR17_valid_signal_gain", Math.Abs(clean - 9.0 / 2020) <= 2e-16, new[] { clean }, "Exact no-error reference.");
        Case(rows, "FR18_error_salience_counterexample", noisy > clean && Math.Abs(noisy - 13.0 / 945) <= 2e-16,
            new[] { noisy }, "Error can raise salience in the presence of a valid signal; no general noise-monotonicity theorem.");
        Case(rows, "FR19_pure_error_erosion", Math.Abs(pure + 1.0 / 500) <= 2e-16, new[] { pure },
            "With zero valid signal, unrecognized error cannot build capacity.");
        var betaSix = Enumerable.Repeat(1.0 / 6, 6).ToArray();
        var betaSeven = new[] { 1.0 / 6, 1.0 / 12, 1.0 / 12, 1.0 / 6, 1.0 / 6, 1.0 / 6, 1.0 / 6 };
        var qvSix = new double[6]; qvSix[0] = 0.5;
        var qvSeven = new double[7]; qvSeven[0] = 0.5;
        double six = CriticalLearning(0.5, Dot(betaSix, qvSix), Dot(betaSix, Enumerable.Repeat(0.01, 6).ToArray())).value.Value;
        double seven = CriticalLearning(0.5, Dot(betaSeven, qvSeven), Dot(betaSeven, Enumerable.Repeat(0.01, 7).ToArray())).value.Value;
        Case(rows, "FR20_taxonomy_counterexample", Near(six, 0.24) && Near(six, seven), new[] { six, seven },
            "Splitting a zero-relevance capacity with split weights need not change the learning threshold.");
        var noEvidence = Period(ExampleInp(DB, new ExParams { sig = 0, noise = 0, chi = 0, lockR = 0 }, stock,
            new[] { 0.55, 0.60 }, 0), false);
        Case(rows, "FR21_zero_signal_validity", noEvidence.v == 0 && noEvidence.ell == 0
            && noEvidence.Cn.All(x => Math.Abs(x - 0.3 * (1 - noEvidence.deltaT)) <= 2e-16),
            new[] { noEvidence.v, noEvidence.ell, noEvidence.Cn[0] }, "Noise-free zero evidence has validity zero and no learning.");
        return rows.ToArray();
    }

    public static string RunInterfaceOnly()
    {
        RunInterfaceChecks();
        return Out.ToString();
    }

    public static string Run(int trials)
    {
        RunRegime(RA, trials);
        RunRegime(RB, trials);
        RunStructure();
        RunDynamics();
        RunConvexity();
        RunSelection(20000);
        RunEdgeCases();
        RunExample();
        RunSensitivity();
        RunGlobal(2000);
        RunHorizon();
        RunInterfaceChecks();
        return Out.ToString();
    }
}
'@

if (-not $DefinitionsOnly -and -not $ScenarioChecksOnly -and -not $SafeguardChecksOnly -and -not $BoundaryChecksOnly) {
    "Environment: Windows PowerShell $($PSVersionTable.PSVersion); .NET CLR $([Environment]::Version); $([Environment]::OSVersion.VersionString); run $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
    "NATIVE_SOURCE_SHA256|$((Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLowerInvariant())"
}
Add-Type -TypeDefinition $src -Language CSharp
if ($DefinitionsOnly) { return }
if ($BoundaryChecksOnly) {
    $checks = @([AfVerify]::BoundaryCases())
    $checks | ConvertTo-Json -Depth 5
    $failed = @($checks | Where-Object { -not $_.Passed })
    if ($checks.Count -ne 21 -or $failed.Count -ne 0) { throw "Boundary checks failed: $($failed.Id -join ', ')" }
    exit 0
}
if ($SafeguardChecksOnly) {
    $checks = @([AfVerify]::SafeguardCases())
    $checks | ConvertTo-Json -Depth 5
    $failed = @($checks | Where-Object { -not $_.Passed })
    if ($checks.Count -eq 0 -or $failed.Count -ne 0) { throw "Safeguard checks failed: $($failed.Id -join ', ')" }
    exit 0
}
if ($ScenarioChecksOnly) {
    $checks = @([AfVerify]::ScenarioCases())
    $checks | ConvertTo-Json -Depth 5
    $failed = @($checks | Where-Object { -not $_.Passed })
    if ($checks.Count -eq 0 -or $failed.Count -ne 0) {
        throw "Scenario checks failed: $($failed.Id -join ', ')"
    }
    exit 0
}
if ($InterfaceChecksOnly) {
    [AfVerify]::RunInterfaceOnly()
    if ([AfVerify]::TestFailures -ne 0) { throw "Interface checks failed: $([AfVerify]::TestFailures)" }
    "INTERFACE CHECKS RESULT: PASS"
    exit 0
}
[AfVerify]::Run(100000)
$scenarioChecks = @([AfVerify]::ScenarioCases())
$scenarioFailures = @($scenarioChecks | Where-Object { -not $_.Passed })
if ($scenarioChecks.Count -eq 0 -or $scenarioFailures.Count -ne 0) {
    throw "Scenario checks failed: $($scenarioFailures.Id -join ', ')"
}
"Scenario checks passed: $($scenarioChecks.Count)"
$safeguardChecks = @([AfVerify]::SafeguardCases())
$safeguardFailures = @($safeguardChecks | Where-Object { -not $_.Passed })
if ($safeguardChecks.Count -eq 0 -or $safeguardFailures.Count -ne 0) { throw "Safeguard checks failed: $($safeguardFailures.Id -join ', ')" }
"Safeguard checks passed: $($safeguardChecks.Count)"
$boundaryChecks = @([AfVerify]::BoundaryCases())
$boundaryFailures = @($boundaryChecks | Where-Object { -not $_.Passed })
if ($boundaryChecks.Count -ne 21 -or $boundaryFailures.Count -ne 0) { throw "Boundary checks failed: $($boundaryFailures.Id -join ', ')" }
"Boundary checks passed: $($boundaryChecks.Count)"

# Values reported in Supplementary Material 1
$expected = @(

    'A|classes|93770,4109,1431,0,684,6',
    'A|status(gaining/stable/eroding by class)|5248/0/88522,150/0/3959,91/0/1340,0/0/0,28/0/656,0/0/6',
    'B|classes|58535,868,30477,1711,4714,3695',
    'B|status(gaining/stable/eroding by class)|15552/0/42983,77/0/791,7554/0/22923,1711/0/0,249/0/4465,60/0/3635',
    'C|N1|N=2000|count=0',
    'C|O3|N=2000|count=5',
    'C|worst|2.22E-015',
    'C|feedback slow/multi|7/5',
    'C|bistable example|theta=0.235,indices=0.273/0.864',
    'S|O1|N=20000|count=11669',
    'S|O2|N=20000|count=1319',
    'W|A|MQ,ell,ellR,v,crit|0.0525,0.0085,0.0012,0.8000,0.0565',
    'W|A|dC1,dCd1|-0.0077,0.0001',
    'W|A|C5,Ccf5,SC5,social5|0.3630,0.3625,0.5519,0.3597',
    'W|A|Ccf5noNoise,maxPeriodGain|0.3626,0.0001',
    'W|A|classes,status|UUXXX,-----',
    'W|A|Csteady|0.0063',
    'W|B|MQ,ell,ellR,v,crit|0.4536,0.3305,0.1851,0.8000,0.0629',
    'W|B|dC1,dCd1|0.0170,0.0251',
    'W|B|C5,Ccf5,SC5,social5|0.4811,0.3612,0.6122,0.6610',
    'W|B|Ccf5noNoise,maxPeriodGain|0.3704,0.0251',
    'W|B|classes,status|GGGGG,+++++',
    'W|B|Csteady|0.7662',
    'W|jensen|h=1,noise=0.0|ystar=0.0000,gap(0.12+/-0.08)=-6.86E-004,gap(0.80+/-0.40)=-1.47E-003',
    'W|jensen|h=1,noise=0.2|ystar=0.0000,gap(0.12+/-0.08)=-3.75E-004,gap(0.80+/-0.40)=-1.51E-003',
    'W|jensen|h=2,noise=0.0|ystar=0.2309,gap(0.12+/-0.08)=8.24E-004,gap(0.80+/-0.40)=-3.87E-003',
    'W|jensen|h=2,noise=0.2|ystar=0.2309,gap(0.12+/-0.08)=9.24E-005,gap(0.80+/-0.40)=-2.96E-003',
    'Z|baseline|A:C5=0.363,cls=UUXXX|B:C5=0.481,gainVsCf=0.120,cls=GGGGG|B-A=0.118',
    'Z|delta0=0.05|A:C5=0.314,cls=UUXXX|B:C5=0.439,gainVsCf=0.117,cls=GGGGG|B-A=0.125',
    'Z|delta0=0.08|A:C5=0.271,cls=UUXXX|B:C5=0.400,gainVsCf=0.113,cls=RRGGG|B-A=0.130',
    'Z|eta=0.3|A:C5=0.363,cls=UUXXX|B:C5=0.421,gainVsCf=0.057,cls=GGGGG|B-A=0.058',
    'Z|eta=0.9|A:C5=0.363,cls=UUXXX|B:C5=0.539,gainVsCf=0.180,cls=GGGGG|B-A=0.176',
    'Z|h=1|A:C5=0.363,cls=UUXXX|B:C5=0.461,gainVsCf=0.106,cls=GGGGG|B-A=0.098',
    'Z|h=3|A:C5=0.363,cls=UUXXX|B:C5=0.494,gainVsCf=0.129,cls=GGGGG|B-A=0.131',
    'Z|theta=1.0|A:C5=0.365,cls=UUXXX|B:C5=0.578,gainVsCf=0.220,cls=GGGGG|B-A=0.214',
    'Z|theta=0.0|A:C5=0.361,cls=UUXXX|B:C5=0.398,gainVsCf=0.036,cls=RRRRR|B-A=0.036',
    'Z|phi=0.0|A:C5=0.360,cls=UUXXX|B:C5=0.470,gainVsCf=0.118,cls=GGGGG|B-A=0.109',
    'Z|phi=1.0|A:C5=0.366,cls=UUXXX|B:C5=0.493,gainVsCf=0.121,cls=GGGGG|B-A=0.127',
    'Z|q=0.25|A:C5=0.363,cls=UUXXX|B:C5=0.421,gainVsCf=0.057,cls=GGGGG|B-A=0.058',
    'Z|q=1.00|A:C5=0.364,cls=UUXXX|B:C5=0.594,gainVsCf=0.238,cls=GGGGG|B-A=0.231',
    'Z|noise=0.05|A:C5=0.363,cls=UUXXX|B:C5=0.505,gainVsCf=0.135,cls=GGGGG|B-A=0.142',
    'Z|noise=0.60|A:C5=0.363,cls=UUXXX|B:C5=0.445,gainVsCf=0.104,cls=GGGGG|B-A=0.083',
    'Z|chi=0.10|A:C5=0.361,cls=UUXXX|B:C5=0.440,gainVsCf=0.095,cls=GGGGG|B-A=0.079',
    'Z|period=5y|A:C5=0.363,cls=X|B:C5=0.401,gainVsCf=0.033,cls=X|B-A=0.039',
    'G|shares(%)|BaboveA=96.2,BaboveCounterfactual=97.1,BcapabilityGainEveryYear=17.9,BcapabilityGainAnyYear=22.4,BlevelRises=37.9,AanyFavourable=0.0',
    'G|noiseFreeCounterfactual|BaboveNoiseFreeCf(%)=87.4,maxPeriodGainB=0.1584',
    'G|spearman(levelChange)|theta=0.49,q=0.43,delta0=-0.43,c0=-0.35,eta=0.30,phi=0.21,disturbance=0.21,chi=-0.16,noise=-0.13,sigma0=-0.11,h=0.06,kappa=0.04,lock=0.00',
    'G|spearman(gainVsCounterfactual)|q=0.58,theta=0.56,eta=0.41,disturbance=0.27,sigma0=-0.21,chi=-0.11,h=0.06,phi=0.06,c0=-0.05,delta0=-0.05,noise=-0.02,kappa=0.02,lock=-0.01',
    'H|A|C0,C30,Ccf30,levelChange,retainedGain,gainPeriods|0.4000,0.2217,0.2212,-1.78E-001,0.0005,0',
    'H|A|classes|UUUUUUUUUUUXXXXXXXXXXXXXXXXXXX',
    'H|A|judgement|no retained-gain verdict (initial/review protections, monitoring, or fixed-baseline performance failed)',
    'H|B|C0,C30,Ccf30,levelChange,retainedGain,gainPeriods|0.4000,0.3373,0.2154,-6.27E-002,0.1218,4',
    'H|B|classes|RRGRRRRGRRRRGRRRXGRXXXXXXXXXXX',
    'H|B|judgement|no retained-gain verdict (initial/review protections, monitoring, or fixed-baseline performance failed)',
    'H|B-mature|C0,C30,Ccf30,levelChange,retainedGain,gainPeriods|0.2761,0.2761,0.1490,0.00E+000,0.1271,0',
    'H|B-mature|classes|XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
    'H|B-mature|judgement|no retained-gain verdict (initial/review protections, monitoring, or fixed-baseline performance failed)',
    'H|mature|dC,dCd,class|0.00E+000,0.0151,CapabilityGain',
    'R|EC0|latentClass,CGI,capacityEvidenceStatus|CapabilityGain,0.0000,InsufficientEvidence',
    'R|cap5|A|d,social,harm,C,latentClass|1.5000,0.4773,0.4500,0.3627,Unacceptable',
    'R|cap5|B|d,social,harm,C,latentClass|1.5000,0.5345,0.4500,0.4013,Unacceptable',
    'R|uncap5|A|d,social,harm,C,latentClass|4.0000,0.3525,1.2000,0.3627,Unacceptable',
    'R|uncap5|B|d,social,harm,C,latentClass|4.0000,0.4099,1.2000,0.4066,Unacceptable',
    'R|compoundDirection|jensenGap|-5.2959097384E-006'

)
# Number of cases on which each randomized test was evaluated
$expectedCounts = @(

    'A|T1|N=100000',
    'A|T2|N=100000',
    'A|T3|N=100000',
    'A|T11|N=100000',
    'A|T12|N=100000',
    'A|T13|N=100000',
    'A|T14|N=100000',
    'A|T15|N=100000',
    'A|T4|N=49984',
    'A|T5|N=49984',
    'A|T6|N=49984',
    'A|T8|N=49984',
    'A|T10|N=49984',
    'A|T9|N=4965',
    'A|T7|N=9482',
    'B|T1|N=100000',
    'B|T2|N=100000',
    'B|T3|N=100000',
    'B|T11|N=100000',
    'B|T12|N=100000',
    'B|T13|N=100000',
    'B|T14|N=100000',
    'B|T15|N=100000',
    'B|T4|N=49847',
    'B|T5|N=49847',
    'B|T6|N=49847',
    'B|T7|N=9509',
    'B|T8|N=49847',
    'B|T10|N=49847',
    'B|T9|N=5081',
    'B|T16|N=1711',
    'P|T17|N=20000',
    'P|T18|N=20000',
    'C|T19|N=10000',
    'C|T20|N=2000',
    'V|T21|N=6968',
    'V|T22|N=10000',
    'S|T23|N=20000',
    'S|T24|N=20000'

)

$actual = [AfVerify]::Reported
$missing = @($expected | Where-Object { $actual -notcontains $_ })
$actualCounts = [AfVerify]::Counts
$countMismatch = @($expectedCounts | Where-Object { $actualCounts -notcontains $_ }) + @($actualCounts | Where-Object { $expectedCounts -notcontains $_ })
$failures = [AfVerify]::TestFailures
""
"== SUMMARY"
"Test failures (Property, Implementation, Consistency, Edge): $failures"
if ($missing.Count -eq 0) { "Reported values reproduced: all $($expected.Count)" }
else { "Reported values NOT reproduced: $($missing.Count) of $($expected.Count)"; $missing | ForEach-Object { "  expected: $_" } }
if ($countMismatch.Count -eq 0) { "Test case counts as reported: all $($expectedCounts.Count)" }
else { "Test case counts NOT as reported:"; $countMismatch | ForEach-Object { "  mismatch: $_" } }
if ($failures -eq 0 -and $missing.Count -eq 0 -and $countMismatch.Count -eq 0) { "RESULT: PASS"; exit 0 } else { "RESULT: FAIL"; exit 1 }
