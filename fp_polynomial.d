//module fp_polynomial;
//
//import std.stdio;
//import std.functional : memoize;
//import std.algorithm : reduce, max, min;
//alias reduce!max maximum;
//
//import fint : F = F16;
//alias F[] FP;
//alias ushort uinteger;
//const N = 16;
//
//version = demo;
//version = naive;
//
//alias size_t[] BP;
//const B = size_t.sizeof << 3;
//static if (B == 32) const BB = 5;
//else static if (B == 64) const BB = 6;
//else static assert (false);
//
//
//FP BP_to_FP(BP f)
//{
//    FP ret;
//    ret.length = f.length << BB;
//    size_t d;
//    foreach (c; f)
//    {
//        foreach (i; 0..B)
//        {
//            if (c >> i & 1)
//            {
//                ret[d] = F(1);
//            }
//            d += 1;
//        }
//    }
//    return ret;
//}
//
//size_t exponent_of_two(size_t n)
//{
//    assert (n);
//    if (n & 1)
//    {
//        return 0;
//    }
//    return (n >> 1).exponent_of_two() + 1;
//}
//
//size_t factorial_exponent_of_two(size_t n)
//{
//    size_t count;
//    foreach (i; 0..n)
//    {
//        count += (i + 1).exponent_of_two();
//    }
//    return count;
//}
//
//bool choose(uint n, uint r)
//{
//    return n.factorial_exponent_of_two == r.factorial_exponent_of_two + (n - r).factorial_exponent_of_two;
//}
//
//size_t[] minimal_polynomial_(uint k)
//{
//    size_t[] ret;
//    foreach (i; 0..k)
//    {
//        if (k.choose(i))
//            ret ~= 1 << i;
//    }
//    return ret ~ (1 << k);
//}
//
//F[] SFT(FP f, uint k)
//{
//    F[] ret;
//    F cur, p, x;
//    ret.length = 1 << k;
//    foreach (i; 0..(1 << k))
//    {
//        cur = F(0);
//        p = F(1);
//        x = F(cast(uinteger)i);
//        foreach (c; f)
//        {
//            cur += c * p;
//            p *= x;
//        }
//        ret[i] = cur;
//    }
//    return ret;
//}
//
//alias memoize!minimal_polynomial_ minimal_polynomial;
//
//FP construct_quotient(size_t[] degree, F[] coefficient)
//{
//    if (degree.length == 0)
//    {
//        return [];
//    }
//    FP ret;
//    ret.length = maximum(degree) + 1;
//    foreach (i, d; degree)
//    {
//        ret[d] = coefficient[i];
//    }
//    return ret;
//}
//
//FP[2] divide_by_sparse(FP f, size_t[] sparse, F const_term)
//{
//    auto divisor_deg = sparse[$ - 1];
//    FP remainder = f.dup;
//    size_t[] qd;
//    F[] qc; 
//    size_t cqd;
//    F cqc;
//    foreach_reverse (i; divisor_deg..remainder.length)
//    {
//        if (!remainder[i])
//        {
//            continue;
//        }
//        cqd = i - divisor_deg;
//        cqc = remainder[i];
//        remainder[cqd] -= cqc * const_term;
//        foreach (j; sparse)
//        {
//            remainder[j + cqd] -= cqc;
//        }
//        qd ~= cqd;
//        qc ~= cqc;
//    }
//    return [construct_quotient(qd, qc), remainder[0..min($, divisor_deg)]];
//}
//
//FP add(FP one, FP another)
//{
//    if (one.length < another.length)
//    {
//        return another.add(one);
//    }
//    FP ret = one.dup;
//    foreach (i, c; another)
//    {
//        ret[i] += c;
//    }
//    return ret;
//}
//
//struct WZC16
//{
//    F[] basis;
//    this (bool dummy)
//    {
//        this.basis = F.basis();
//    }
//    F pi(uinteger x)
//    {
//        F ret;
//        foreach (i, b; this.basis)
//        {
//            if (x >> i & 1)
//            {
//                ret += b;
//            }
//        }
//        return ret;
//    }
//    /*
//    let p(x) = x^2 + x.
//    then minimal polynomial of F_i = p^i(x).
//    */
//    FP[2] divide_by_minimal(FP f, uint i, uinteger j)
//    {
//        return f.divide_by_sparse(minimal_polynomial(i), this.pi(j));
//    }
//    F[] FFT(FP f)
//    {
//        "transforming %s ...".writefln(f);
//        auto ret = this.FFT_(f, N, 0);
//        "returning %s".writefln(ret);
//        return ret;
//    }
//    private F ev(FP f, uint j)
//    {
//        return f[0];
//    }
//    /*
//input: f, i, j
//output: f[j << i .. (j+1) << i)
//    */
//    F[] FFT_(FP f, uint i, uinteger j)
//    {
//        if (i == 0)
//        {
//            return [this.ev(f, j)];
//        }
//        i -= 1;
//        j <<= 1;
//        auto qr = this.divide_by_minimal(f, i, j);
//        auto f0 = qr[1];
//        auto f1 = qr[0].add(qr[1]);
//        return
//            FFT_(f0, i, j) ~ FFT_(f1, i, j | 1);
//    }
//}
//
//version (demo) unittest
//{
//    FP f;
//    FP g;
//    FP h;
//    foreach (i; 0..10)
//    {
//        f ~= F(cast(uinteger)i);
//        g ~= F(cast(uinteger)(10 - i));
//        h ~= f[$-1] + g[$-1];
//    }
//    auto ffter = WZC16(false);
//    auto Tf = ffter.FFT_(f, 4, 0), Tg = ffter.FFT_(g, 4, 0), Th = ffter.FFT_(h, 4, 0);
//    foreach (i, v; Th)
//    {
//        assert (Tf[i] + Tg[i] == v);
//    }
//    version (naive)
//    {
//        auto s = SFT(f, N);
//        foreach (i; 0..(1 << 4))
//        {
//            assert (Tf[i] == s[ffter.pi(cast(uinteger)i).v]);
//        }
//    }
//}
