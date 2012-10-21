import std.functional : memoize;
import std.algorithm : min, max, reduce;
debug import std.stdio : writefln;
alias reduce!max maximum;

alias size_t[] BP;
const size_t B = size_t.sizeof << 3;
static if (B == 8) const size_t BB = 3;
else static if (B == 16) const size_t BB = 4;
else static if (B == 32) const size_t BB = 5;
else static if (B == 64) const size_t BB = 6;
else static if (B == 128) const size_t BB = 7;
else static assert (false);
const size_t mask = B - 1;

/** number of trailing zeros of n */
size_t exponent_of_two(size_t n)
{
    assert (n);
    size_t ret;
    while ((n & 1) == 0)
    {
        n >>= 1;
        ret += 1;
    }
    return ret;
}

/** number of trailing zeros of n! */
size_t factorial_exponent_of_two(size_t n)
{
    size_t count;
    foreach (i; 0..n)
    {
        count += (i + 1).exponent_of_two();
    }
    return count;
}

/** binomial coefficient modulo 2 */
bool choose(size_t n, size_t r)
{
    return n.factorial_exponent_of_two == r.factorial_exponent_of_two + (n - r).factorial_exponent_of_two;
}

/** minimal polynomial of F_k */
size_t[] minimal_polynomial_(uint k)
{
    size_t[] ret;
    foreach (i; 0..k)
        if (k.choose(i))
            ret ~= 1 << i;
    return ret ~ (1 << k);
}
alias memoize!minimal_polynomial_ minimal_polynomial; /// ditto

/** minimize the length of array according to the polynomial */
BP simp(BP f)
{
    foreach_reverse (i, c; f)
    {
        if (c)
        {
            return f[0..(i+1)];
        }
    }
    return [];
}

/** tool for degree */
ptrdiff_t highest_bit(size_t x)
{
    ptrdiff_t ret = -1;
    while (x)
    {
        ret += 1;
        x >>= 1;
    }
    return ret;
}

/** degree of binary polynomial */
ptrdiff_t degree(BP f)
{
    foreach_reverse (i, c; f)
    {
        if (c)
        {
            return (i << BB) + c.highest_bit();
        }
    }
    return -1;
}

/** B[x] -> F[x] */
F[] to_FP(F)(BP f)
{
    F[] ret;
    ret.length = f.degree + 1;
    auto one = F(1);
    size_t d;
    foreach (c; f)
    {
        foreach (i; 0..B)
        {
            if (c >> i & 1)
            {
                ret[d] = one;
            }
            d += 1;
        }
    }
    return ret;
}


mixin template FFTTOOLS(F, U)
{
    alias F[] FP;
    BP fftprod(BP f, BP g)
    {
        return f.simp().to_FP!F().fftprod(g.simp().to_FP!F()).to_BP.simp();
    }
    FP fftprod(FP f, FP g)
    {
        auto d = f.degree(), e = g.degree();
        if (d < 0 || e < 0)
        {
            return [];
        }
        uint i = 0;
        auto n = d + e;
        while (1 << i <= n)
        {
            i += 1;
        }
        // 1 << i > f.deg + g.deg
        auto Tf = fft.FFT(f, i);
        auto Tg = fft.FFT(g, i);
        //assert (equality!(F, FP)(fft.IFFT(Tf, i), f));
        //assert (equality!(F, FP)(fft.IFFT(Tg, i), g));
        auto Tp = Tf.pointwise_product(Tg);
        return fft.IFFT(Tp, i);
    }

    bool equality(F, FP : F[])(FP x, FP y)
    {
        foreach (i, c; x)
        {
            if (c != y[i])
                return false;
        }
        return true;
    }

    void pad(ref FP f, uint i)
    {
        f.length = 1 << i;
    }

    F[] pointwise_product(F[] x, F[] y) /// pointwise product of two arrays of same length
    {
        F[] ret;
        ret.length = x.length;
        assert (ret.length == y.length);
        foreach (i, c; x)
        {
            ret[i] = c * y[i];
        }
        return ret;
    }

    version(none){
    /** simplify F[x] by omitting zeros */
    FP simp(FP f)
    {
        return f[0..(f.degree() + 1)];
    }}

    /// degree
    ptrdiff_t degree(FP f)
    {
        foreach_reverse (i, c; f)
        {
            if (c)
            {
                return i;
            }
        }
        return -1;
    }

    /** F[x] -> B[x] */
    BP to_BP(FP f)
    {
        BP ret;
        ret.length = (f.length >> BB) + 1;
        size_t bit = 1;
        foreach (i, c; f)
        {
            if (!c)
            {
                continue;
            }
            assert (c == F(1));
            ret[i >> BB] |= bit << (i & mask);
        }
        return ret;
    }

    /** slow fourier transform: O(n << k) */
    F[] SFT(FP f, uint k)
    {
        F[] ret;
        F cur, p, x;
        ret.length = 1 << k;
        foreach (i; 0..(1 << k))
        {
            cur = F(0);
            p = F(1);
            x = F(cast(U)i);
            foreach (c; f)
            {
                cur += c * p;
                p *= x;
            }
            ret[i] = cur;
        }
        return ret;
    }

    version(none){/** [d], [c] -> [d-th element is c] */
    FP construct_quotient(size_t[] degree, F[] coefficient)
    {
        if (degree.length == 0)
        {
            return [];
        }
        FP ret;
        ret.length = maximum(degree) + 1;
        foreach (i, d; degree)
        {
            ret[d] = coefficient[i];
        }
        return ret;
    }}

    /** divide f by ((x^j : j in sparse).sum + const_term) */
    FP[2] divide_by_sparse(FP f, size_t[] sparse, F const_term)
    {
        debug (fft) "dividing %s by sparse%s + %s".writefln(f, sparse, const_term);
        auto divisor_deg = sparse[$ - 1];
        FP remainder = f.dup;
        FP quotient;
        quotient.length = f.length - divisor_deg; // q.d = f.d - dd < f.len - dd
        foreach_reverse (i; divisor_deg..remainder.length)
        {
            if (!remainder[i])
            {
                continue;
            }
            auto cqc = remainder[i];
            size_t cqd = i - divisor_deg;
            quotient[cqd] = cqc;
            remainder[cqd] -= cqc * const_term;
            foreach (j; sparse)
            {
                remainder[j + cqd] -= cqc;
            }
        }
        debug (fft) "returning %s ... %s".writefln(quotient, remainder[0..min($, divisor_deg)]);
        return [quotient, remainder[0..divisor_deg]];
    }

    /** multiply f  by ((x^j : j in sparse).sum + const_term) */
    FP multiply_by_sparse(FP f, size_t[] sparse, F const_term)
    {
        //debug "multiplying %d-element polynomial by sparse%s + %s ...".writefln(f.length, sparse, const_term);
        auto ret = f.dup;
        foreach (i, c; f)
        {
            ret[i] *= const_term;
        }
        ret.length += sparse[$ - 1];
        foreach (i, c; f)
        {
            foreach (j; sparse)
            {
                ret[i + j] += c;
            }
        }
        //debug "returning %d-element polynomial.".writefln(ret.length);
        return ret;
    }

    /** F[x] + F[x] */
    FP add(FP one, FP another)
    {
        assert (another.length <= one.length);
        FP ret = one.dup;
        foreach (i, c; another)
        {
            ret[i] += c;
        }
        return ret;
    }
}


mixin template WZC(F, U, uint n)
{
    static F[] basis = F.cantor_basis(); /** Cantor basis */
    /** pi: (e[i] << i).sum -> (e[i] b[i]).sum */
    F pi(U x)
    {
        debug (fft) "evaluating pi(%d) ".writef(x);
        F ret;
        foreach (i, b; this.basis)
        {
            if (x >> i & 1)
            {
                ret += b;
            }
        }
        debug (fft) "= %s".writefln(ret);
        return ret;
    }
    /** divide f by s[i]-pi(j) */
    FP[2] divide_by_minimal(FP f, uint i, U j)
    {
        return f.divide_by_sparse(minimal_polynomial(i), pi(j));
    }
    /** perform FFT from top level */
    F[] FFT(FP f, uint i = n)
    {
        pad(f, i);
        debug (fft) "transforming %s ...".writefln(f);
        auto ret = this.FFT_(f, i);
        debug (fft) "returning %s".writefln(ret);
        return ret;
    }
    /** perform FFT recursively.
    Params:
    f, i, j
    Returns:
    [j << i .. (j+1) << i) . pi . f
    */
    F[] FFT_(FP f, uint i, U j = 0)
    {
        assert (f.length == 1 << i);
        debug (fft) "%s[%d..%d)".writefln(f, j << i, (j+1) << i);
        if (i == 0)
        {
            return [f[0]];
        }
        i -= 1;
        j <<= 1;
        auto qr = this.divide_by_minimal(f, i, j);
        auto ret = FFT_(qr[1], i, j) ~ FFT_(qr[0].add(qr[1]), i, j | 1);
        debug (fft) "returning %s.".writefln(ret);
        return ret;
    }

    FP mutiply_by_minimal(FP f, uint i, U j)
    {
        return f.multiply_by_sparse(minimal_polynomial(i), pi(j));
    }
    FP IFFT(F[] f, uint i = n)
    {
        debug (fft) "inverse-transforming %s ...".writefln(f);
        auto ret = this.IFFT_(f, i);
        debug (fft) "returning %s".writefln(ret);
        return ret;
    }
    /** Perform IFFT recursively.
    Params:
    values = [j << i .. (j+1) << i) . pi . f
    Returns:
    f mod s[i+1] - pi(j)
    */
    FP IFFT_(F[] values, uint i, U j = 0)
    {
        assert (values.length == 1 << i);
        if (i == 0)
        {
            return [values[0]];
        }
        i -= 1;
        j <<= 1;
        auto r = IFFT_(values[0..($ >> 1)], i, j);
        //debug "r.length = %d".writefln(r.length);
        auto q = IFFT_(values[($ >> 1)..$], i, j | 1).add(r);
        //debug "r.length = %d".writefln(r.length);
        return this.mutiply_by_minimal(q, i, j).add(r);
    }
}
