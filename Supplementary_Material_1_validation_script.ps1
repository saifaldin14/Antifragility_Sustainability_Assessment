# Verification script for Supplementary Material 1:
# "A Dynamic Mathematical Model for Antifragile Sustainability Assessment".
# Run with: powershell -NoProfile -ExecutionPolicy Bypass -File Supplementary_Material_1_validation_script.ps1
$ErrorActionPreference = 'Stop'
$src = @'
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

public static class AfVal
{
    // Model dimensions used in the randomized test
    const int J = 15, K = 6, m = 3, n = 4, p = 4, r = 3, NG = 3;
    const int iO = 2, iR = 4;                      // positions of O and R in the capacity vector
    const double EpsS = 0.01, EpsC = 0.0, Eps0 = 0.05;
    const double TolBound = 1e-12, TolZero = 1e-15, TolConv = 1e-6;

    static Random Rng;
    static List<string> order;
    static Dictionary<string, long[]> res;
    static Dictionary<string, string> kind;
    static StringBuilder Out = new StringBuilder();

    static double U(double a, double b) { return a + (b - a) * Rng.NextDouble(); }
    static double[] Vec(int k, double a, double b) { var v = new double[k]; for (int i = 0; i < k; i++) v[i] = U(a, b); return v; }
    static double[] Simplex(int k) { var v = new double[k]; double s = 0; for (int i = 0; i < k; i++) { v[i] = -Math.Log(1 - Rng.NextDouble()); s += v[i]; } for (int i = 0; i < k; i++) v[i] /= s; return v; }
    static double[,] Mat(int a, int b, double lo, double hi) { var x = new double[a, b]; for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) x[i, j] = U(lo, hi); return x; }
    static double[] Mul(double[,] A, double[] v) { int a = A.GetLength(0), b = A.GetLength(1); var o = new double[a]; for (int i = 0; i < a; i++) { double s = 0; for (int j = 0; j < b; j++) s += A[i, j] * v[j]; o[i] = s; } return o; }
    static double Dot(double[] a, double[] b) { double s = 0; for (int i = 0; i < a.Length; i++) s += a[i] * b[i]; return s; }
    static double Norm(double[] v) { return Math.Sqrt(Dot(v, v)); }
    static double Clip(double x) { return Math.Max(0, Math.Min(1, x)); }

    // Model components (Equations S13-S14, S17-S18, S23, S32-S33, S37)
    static double Ell(double MQ, double rho, double gamma, double ne, double s0) { return MQ * rho * gamma * ne / (ne + s0); }
    static double[] Response(double ell, double[,] A, double[] e, double kappa, double O, double R)
    {
        double ne = Norm(e); int rr = A.GetLength(0), pp = A.GetLength(1);
        var eh = new double[pp]; if (ne > 0) for (int i = 0; i < pp; i++) eh[i] = e[i] / ne;
        var us = Mul(A, eh); for (int i = 0; i < rr; i++) us[i] *= ell;
        double cap = kappa * (O + R), nu = Norm(us);
        if (nu > cap && nu > 0) for (int i = 0; i < rr; i++) us[i] *= cap / nu;
        return us;
    }
    static double[] NextC(double[] C, double eta, double ell, double[] q, double delta)
    {
        var o = new double[C.Length];
        for (int k = 0; k < C.Length; k++) o[k] = C[k] + eta * ell * q[k] * (1 - C[k]) - delta * C[k];
        return o;
    }
    static double Afi(double dC, double MQ, double EQ, bool SB, double nd)
    {
        double psi = dC > 0 ? MQ * EQ * (SB ? 1 : 0) : 1;
        return psi * dC / (nd + Eps0);
    }
    static int ClassOf(bool SB, double SC, double SCn, double SCmin, double dC)
    {
        if (!SB) return 0;
        if (SCn < SC - EpsS) return 1;
        if (dC > EpsC) return 2;
        return SCmin >= SC - EpsS ? 3 : 4;
    }

    static int SelectBest(double[] V, bool[] feas)
    {
        int best = -1; for (int i = 0; i < V.Length; i++) if (feas[i] && (best < 0 || V[i] > V[best])) best = i;
        return best;
    }

    static void Check(string name, string k, bool ok)
    {
        if (!res.ContainsKey(name)) { res[name] = new long[2]; order.Add(name); kind[name] = k; }
        res[name][0]++; if (!ok) res[name][1]++;
    }
    static void Reset() { order = new List<string>(); res = new Dictionary<string, long[]>(); kind = new Dictionary<string, string>(); }
    static string Ub(long nObs, long fail)
    {
        if (fail > 0) return "n/a";
        double ub = 1 - Math.Pow(0.05, 1.0 / nObs);
        return ub.ToString("E2");
    }
    static void Dump(string header)
    {
        Out.AppendLine("== " + header);
        foreach (var nm in order)
            Out.AppendLine(string.Format("{0}|{1}|N={2}|count={3}|rate={4:P2}|UB95={5}", kind[nm], nm, res[nm][0], res[nm][1], (double)res[nm][1] / res[nm][0], Ub(res[nm][0], res[nm][1])));
    }

    class Regime
    {
        public string Name; public int Seed;
        public double SLo, TauHi, EtaLo, DeltaCap, MQLo, RGLo, BLo, BHi, GLo, GHi, HLo, HHi;
    }

    static void RunRegime(Regime R, int trials)
    {
        Rng = new Random(R.Seed); Reset();
        long[] cls = new long[5];
        for (int t = 0; t < trials; t++)
        {
            var S = Vec(J, R.SLo, 1); var w = Simplex(J); var tau = Vec(J, 0, R.TauHi);
            var crit = new bool[J]; for (int j = 0; j < J; j++) crit[j] = Rng.NextDouble() < 0.3;
            var C = Vec(K, 0, 1); var beta = Simplex(K); var q = Simplex(K);
            double eta = U(R.EtaLo, 1), delta = U(0, Math.Min(R.DeltaCap, 1 - eta));
            var d = Vec(m, 0, 1.5); var E = Mat(n, m, 0, 1); var x = Mul(E, d);
            var M = Mat(p, n, 0, 1.5); var e = Mul(M, x); var noise = Vec(p, -0.2, 0.2); for (int i = 0; i < p; i++) e[i] += noise[i];
            double MQ = U(R.MQLo, 1) * U(R.MQLo, 1) * U(R.MQLo, 1) * U(R.MQLo, 1);
            double rho = Rng.NextDouble() < 0.1 ? 0 : U(R.RGLo, 1);
            double gamma = Rng.NextDouble() < 0.1 ? 0 : U(R.RGLo, 1);
            double s0 = U(0.1, 2), ne = Norm(e);

            double ell = Ell(MQ, rho, gamma, ne, s0);
            Check("T1 Effective learning lies in [0,1]", "Property", ell >= 0 && ell <= 1);

            var A = Mat(r, p, -1, 1); double kappa = U(0, 1);
            var u = Response(ell, A, e, kappa, C[iO], C[iR]);
            Check("T2 Implemented response respects kappa(O+R)", "Implementation", Norm(u) <= kappa * (C[iO] + C[iR]) + TolBound);

            double bHi = U(R.BLo, R.BHi);
            var B = Mat(J, n, 0, bHi); var G = Mat(J, r, R.GLo, R.GHi); var Om = Vec(J, -0.05, 0.05);
            var Bx = Mul(B, x); var Gu = Mul(G, u);
            var Smin = new double[J]; var Sn = new double[J]; bool sOk = true;
            for (int j = 0; j < J; j++) { Smin[j] = Clip(S[j] - Bx[j]); Sn[j] = Clip(S[j] - Bx[j] + Gu[j] + Om[j]); if (Sn[j] < 0 || Sn[j] > 1) sOk = false; }
            Check("T3 Performance remains in [0,1]", "Implementation", sOk);

            var Cn = NextC(C, eta, ell, q, delta);
            bool cOk = true, cTight = true;
            for (int k = 0; k < K; k++)
            {
                double a = eta * ell * q[k];
                if (Cn[k] < -TolBound || Cn[k] > 1 + TolBound) cOk = false;
                if (Cn[k] < a - TolBound || Cn[k] > 1 - delta + TolBound) cTight = false;
            }
            Check("T4 Capacity remains in [0,1] (Property 1)", "Property", cOk);
            Check("T5 Capacity lies in [a_k, 1-delta] (Property 1, tight bound)", "Property", cTight);
            double dC = Dot(beta, Cn) - Dot(beta, C);

            var qc = new double[K]; for (int k = 0; k < K; k++) qc[k] = q[k] * (1 - C[k]);
            double lhs = eta * ell * Dot(beta, qc), rhs = delta * Dot(beta, C);
            if (Math.Abs(lhs - rhs) > TolBound) Check("T6 Learning-decay condition predicts sign of capability change (Property 3)", "Property", (dC > 0) == (lhs > rhs));

            if (rho == 0 || gamma == 0) Check("T7 No capability gain without retention or connectivity (Property 2)", "Property", ell == 0 && dC <= TolZero);

            double MQ2 = MQ + U(0, 1 - MQ), rho2 = rho + U(0, 1 - rho), gamma2 = gamma + U(0, 1 - gamma), scale = U(1, 3);
            double ell2 = Ell(MQ2, rho2, gamma2, ne * scale, s0);
            Check("T8 Learning is non-decreasing in its inputs", "Property", ell2 >= ell - TolZero);
            double dC2 = Dot(beta, NextC(C, eta, ell2, q, delta)) - Dot(beta, C);
            Check("T9 Capability change is non-decreasing in learning", "Property", dC2 >= dC - TolBound);

            var pi = Vec(NG, 0.5, 2); var Hm = Mat(NG, n, 0, 0.5); var h = Mul(Hm, x);
            double H = Dot(pi, h), Hmax = U(R.HLo, R.HHi) * H;
            bool SB = H <= Hmax; for (int j = 0; j < J; j++) if (crit[j] && Sn[j] < tau[j]) SB = false;

            double SC = Dot(w, S), SCn = Dot(w, Sn), SCmin = Dot(w, Smin);
            bool notDown = SCn >= SC - EpsS;
            bool[] k5 = {
                !SB,
                SB && !notDown,
                SB && notDown && dC > EpsC,
                SB && notDown && dC <= EpsC && SCmin >= SC - EpsS,
                SB && notDown && dC <= EpsC && SCmin < SC - EpsS };
            int nTrue = k5.Count(z => z);
            Check("T10 Exactly one response class is assigned (Property 6)", "Property", nTrue == 1);
            for (int i = 0; i < 5; i++) if (k5[i]) cls[i]++;
            bool strong = dC > EpsC && SCn >= SC - EpsS && SB;
            Check("T11 Strong condition (S28) coincides with the antifragile class", "Consistency", strong == k5[2]);

            double EQ = U(0, 1), nd = Norm(d);
            double afi = Afi(dC, MQ, EQ, SB, nd);
            Check("T12 Index sign is consistent with capability change (Property 7)", "Property", ((afi > 0) == (dC > 0 && SB && MQ * EQ > 0)) && ((afi < 0) == (dC < 0)));
            Check("T13 Index is bounded by 1/eps0", "Property", Math.Abs(afi) <= 1 / Eps0 + TolBound);
            if (k5[2]) Check("T14 Antifragile class implies a positive index when MQ*EQ > 0", "Consistency", MQ * EQ <= 0 || afi > 0);

            // Preliminary formulation, evaluated on the same scenario (differential test)
            double ellOld = MQ * rho * gamma * ne; bool over = false;
            for (int k = 0; k < K; k++) if (C[k] + eta * rho * gamma * ellOld * q[k] - delta * C[k] > 1) over = true;
            Check("D1 Preliminary: a capacity exceeds 1", "Defect", !over);
            double dSC = SCn - SC;
            bool oF = dSC < 0 && dC <= 0, oRb = Math.Abs(dSC) <= EpsS && Math.Abs(dC) <= 1e-3, oRs = Math.Abs(dSC) <= EpsS && dC >= 0, oA = dC > 0 && SB;
            int oN = (oF ? 1 : 0) + (oRb ? 1 : 0) + (oRs ? 1 : 0) + (oA ? 1 : 0);
            Check("D2 Preliminary: outcome assigned no class", "Defect", oN != 0);
            Check("D3 Preliminary: outcome assigned more than one class", "Defect", oN <= 1);
            if (dC < 0)
            {
                double full = dC / (nd + Eps0), afiOld = full * MQ * EQ * (SB ? 1 : 0);
                Check("D4 Preliminary: capability loss reported as zero", "Defect", afiOld != 0);
                Check("D5 Preliminary: capability loss understated", "Defect", afiOld == 0 || Math.Abs(afiOld) >= Math.Abs(full));
            }
        }
        Dump("REGIME " + R.Name + " (seed " + R.Seed + ", " + trials + " scenarios)");
        string[] names = { "Unacceptable", "Fragile", "Antifragile", "Robust", "Resilient" };
        Out.AppendLine("classes|" + string.Join("|", Enumerable.Range(0, 5).Select(i => names[i] + "=" + cls[i] + " (" + (100.0 * cls[i] / trials).ToString("F2") + "%)")));
    }

    static void RunConvergence(int trials)
    {
        Rng = new Random(20261008); Reset();
        double worst = 0;
        for (int t = 0; t < trials; t++)
        {
            double eta = U(0, 1), delta = U(0.001, 1 - eta), ell = U(0, 1), qk = U(0, 1), c = U(0, 1);
            double a = eta * ell * qk;
            for (int s = 0; s < 20000; s++) c = c + a * (1 - c) - delta * c;
            double err = Math.Abs(c - a / (a + delta)); worst = Math.Max(worst, err);
            Check("T15 Capacity converges to a/(a+delta) (Property 4)", "Property", err <= TolConv);
        }
        Dump("CONVERGENCE (seed 20261008, " + trials + " trajectories, 20,000 steps each)");
        Out.AppendLine("worst absolute error|" + worst.ToString("E2"));
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
                feas[i] = minCrit >= tauC && Hi <= Hmax && !ir;
            }
            int best = SelectBest(V, feas);
            var ranked = Enumerable.Range(0, NA).Where(i => feas[i]).OrderByDescending(i => V[i]).ToList();
            int indep = ranked.Count > 0 ? ranked[0] : -1;
            Check("T16 Selected alternative is feasible (Property 8)", "Property", best < 0 || feas[best]);
            Check("T17 Selection agrees with independent ranking of feasible alternatives", "Implementation", best == indep);
            int ub = 0; for (int i = 1; i < NA; i++) if (V[i] > V[ub]) ub = i;
            Check("D6 Unconstrained maximum is infeasible", "Defect", feas[ub]);
            Check("O1 No feasible alternative exists", "Observation", best >= 0);
        }
        Dump("ALTERNATIVE SELECTION (seed 20261009, " + trials + " decision problems, 5 alternatives each)");
    }

    static void RunEdgeCases()
    {
        Rng = new Random(20261010); Reset();
        var beta = Simplex(K); var q = Simplex(K); var C = Vec(K, 0, 1); var A = Mat(r, p, -1, 1);
        // E1 zero evidence
        double l1 = Ell(1, 1, 1, 0, 1); var u1 = Response(l1, A, new double[p], 1, 1, 1);
        double dc1 = Dot(beta, NextC(C, 0.5, l1, q, 0.2)) - Dot(beta, C);
        Check("E1 Zero evidence: no learning, no response, capacity decays at rate delta", "Edge", l1 == 0 && Norm(u1) == 0 && Math.Abs(dc1 + 0.2 * Dot(beta, C)) < TolBound);
        Check("E2 Zero retention: no learning", "Edge", Ell(1, 0, 1, 5, 1) == 0);
        Check("E3 Zero connectivity: no learning", "Edge", Ell(1, 1, 0, 5, 1) == 0);
        var ones = Enumerable.Repeat(1.0, K).ToArray(); var zeros = new double[K];
        Check("E4 Capacities at maximum, full learning, no decay: no overshoot", "Edge", NextC(ones, 1, 1, q, 0).All(v => Math.Abs(v - 1) < TolBound));
        Check("E5 Capacities at zero, no learning, no decay: unchanged", "Edge", NextC(zeros, 0.7, 0, q, 0).All(v => v == 0));
        var qk = new double[K]; qk[0] = 1;
        bool e6 = NextC(ones, 0.6, 1, qk, 0.4).Concat(NextC(zeros, 0.6, 1, qk, 0.4)).Concat(NextC(C, 0.6, 1, qk, 0.4)).All(v => v >= -TolBound && v <= 1 + TolBound);
        Check("E6 Boundary parameters eta + delta = 1: capacity remains in [0,1]", "Edge", e6);
        Check("E7 No optionality or reversibility (O = R = 0): no response", "Edge", Norm(Response(1, A, Vec(p, -1, 1), 1, 0, 0)) == 0);
        double a8 = Afi(0.01, 1, 1, true, 0);
        Check("E8 No disturbance (||d|| = 0): index finite and equal to dC/eps0", "Edge", !double.IsInfinity(a8) && Math.Abs(a8 - 0.01 / Eps0) < TolBound);
        var bigE = Vec(p, 1e5, 1e6); double l9 = Ell(1, 1, 1, Norm(bigE), 1);
        var C9 = NextC(C, 0.9, l9, q, 0.1); var S9 = Vec(J, 0, 1).Select(v => Clip(v - 1e6)).ToArray();
        Check("E9 Extreme disturbance: learning, capacity, and performance remain bounded", "Edge", l9 <= 1 && C9.All(v => v >= 0 && v <= 1) && S9.All(v => v == 0));
        Check("E10 Capability gain with safeguard breach: unacceptable class and zero index", "Edge", ClassOf(false, 0.5, 0.6, 0.6, 0.05) == 0 && Afi(0.05, 1, 1, false, 1) == 0);
        Check("E11 Capability loss with safeguard breach: loss reported in full", "Edge", Math.Abs(Afi(-0.05, 1, 1, false, 1) - (-0.05 / (1 + Eps0))) < TolBound);
        Check("E12 Capability loss with no monitoring (MQ = 0): loss reported in full", "Edge", Math.Abs(Afi(-0.05, 0, 1, true, 1) - (-0.05 / (1 + Eps0))) < TolBound);
        Check("E13 Capability gain with no monitoring (MQ = 0): zero index", "Edge", Afi(0.05, 0, 1, true, 1) == 0);
        Check("E14 Robust: no trough, no capability gain", "Edge", ClassOf(true, 0.5, 0.5, 0.5, 0) == 3);
        Check("E15 Resilient: trough then recovery, no capability gain", "Edge", ClassOf(true, 0.5, 0.5, 0.3, 0) == 4);
        Check("E16 Fragile: no recovery, even with capability gain", "Edge", ClassOf(true, 0.5, 0.3, 0.3, 0.05) == 1);
        Check("E17 Antifragile: trough, recovery, and capability gain", "Edge", ClassOf(true, 0.5, 0.5, 0.3, 0.05) == 2);
        Check("E18 Empty feasible set: no alternative selected", "Edge", SelectBest(new double[] { 0.9, 0.5, 0.1 }, new bool[] { false, false, false }) == -1);
        Dump("EDGE CASES (deterministic)");
    }

    public static string Run(int trials)
    {
        var A = new Regime { Name = "A (stress)", Seed = 20261006, SLo = 0, TauHi = 0.3, EtaLo = 0, DeltaCap = 1, MQLo = 0, RGLo = 0, BLo = 0.15, BHi = 0.15, GLo = -0.2, GHi = 0.5, HLo = 0, HHi = 2 };
        var B = new Regime { Name = "B (coverage)", Seed = 20261007, SLo = 0.2, TauHi = 0.15, EtaLo = 0.2, DeltaCap = 0.15, MQLo = 0.6, RGLo = 0.3, BLo = 0.001, BHi = 0.12, GLo = -0.1, GHi = 0.6, HLo = 0.5, HHi = 3 };
        RunRegime(A, trials);
        RunRegime(B, trials);
        RunConvergence(10000);
        RunSelection(20000);
        RunEdgeCases();
        return Out.ToString();
    }
}
'@
Add-Type -TypeDefinition $src -Language CSharp
[AfVal]::Run(200000)
