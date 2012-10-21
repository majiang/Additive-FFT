module metafint;

import std.stdio : write, writeln;
import std.conv : to;
import std.array : replace, replicate;
debug import std.random;


private struct accumulator
{
    string data;
    string type_name;
    uint indented;
    this (string type_name)
    {
        this.type_name = type_name;
    }
    void indent()
    {
        this.append("{");
        indented += 1;
    }
    void outdent()
    {
        indented -= 1;
        this.append("}");
    }
    void append_indent()
    {
        this.raw_append("    ".replicate(indented));
    }
    void raw_append(string x)
    {
        this.data ~= x;
    }
    void append_newline()
    {
        this.raw_append("\n");
    }
    void indented_append(string x)
    {
        this.indented += 1;
        this.append(x);
        this.indented -= 1;
    }
    void append(string x)
    {
        this.append_indent();
        this.raw_append(x.replace("type_name", type_name));
        this.append_newline();
    }
    string wrap_integer (ulong i)
    {
        string ret = i.to!string;
        if (type_name == "ulong") return ret ~ "UL";
        if (type_name == "uint") return ret;
        if (type_name == "ushort") return ret;
        if (type_name == "ubyte") return ret;
        assert (false);
    }
}

private string get_type_name(uint n)
{
    if (n <= 8) return "ubyte";
    if (n <= 16) return "ushort";
    if (n <= 32) return "uint";
    if (n <= 64) return "ulong";
    if (n <= 128) return "ucent";
    assert (false);
}

private void left_shift(ref ulong pp, ulong p, uint n)
{
    if (pp >> (n - 1) & 1)
    {
        pp = (pp << 1) ^ p;
    }
    else
    {
        pp <<= 1;
    }
    if (n != 64)
    {
        pp &= (1UL << n) - 1;
    }
}

private bool is_not_bits_of_type(uint n)
{
    if (n == 8) return false;
    if (n == 16) return false;
    if (n == 32) return false;
    if (n == 64) return false;
    if (n == 128) return false;
    return true;
}

mixin template declare_prod(T, uint n, ulong p)
{
    private static T prod(T a, T b)
    {
        T ret, higher;
        if (a & 1)
        {
            ret = b;
        }
        mixin (gen_prod_internal(n, p));
    }
}

string gen_prod_internal(uint n, ulong p)
{
    accumulator a = accumulator(get_type_name(n));
    a.indented += 1;
    foreach (i; 1..n)
    {
        a.append("if (a & " ~ a.wrap_integer(1UL << i) ~ ")");
        a.indent;{
            a.append("ret ^= b << " ~ i.to!string ~ ";");
            a.append("higher ^= b >> " ~ (n - i).to!string ~ ";");
        }a.outdent;
    }
    if (n.is_not_bits_of_type())
    {
        a.append("ret &= " ~ ((1UL << n) - 1).to!string ~ ";");
    }
    auto pp = p;
    foreach (i; 0..(n - 1))
    {
        a.append("if (higher & " ~ a.wrap_integer(1UL << i) ~ ")");
        a.indented_append("ret ^= " ~ a.wrap_integer(pp) ~ ";");
        pp.left_shift(p, n);
    }
    a.append("return ret;");
    return a.data;
}



string gen_squar(uint n, ulong p)
{
    string type_name = get_type_name(n);
    accumulator a;
    a.type_name = type_name;
    a.append("private static type_name squar(type_name a)");
    a.indent;{
        a.append("type_name ret;");
        ulong pp = 1;
        foreach (i; 0..n)
        {
            a.append("if (a & " ~ a.wrap_integer(1UL << i) ~ ")");
            a.indented_append("ret ^= " ~ a.wrap_integer(pp) ~ ";");
            pp.left_shift(p, n);
            pp.left_shift(p, n);
        }
        a.append("return ret;");
    }a.outdent;
    return a.data;
}

version (generation_large) unittest
{
    gen_prod_internal(64, (1UL << 61) + (1UL << 34) + (1 << 9) + (1 << 0)).writeln();
    gen_squar(64, (1UL << 61) + (1UL << 34) + (1 << 9) + (1 << 0)).writeln();
}
version (generation) unittest
{
    gen_prod_internal(32, (1 << 22) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_squar(32, (1 << 22) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(24, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_squar(24, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(16, (1 << 12) + (1 << 3) + (1 << 1) + (1 << 0)).writeln();
    gen_squar(16, (1 << 12) + (1 << 3) + (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(12, (1 << 10) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_squar(12, (1 << 10) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(8, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_squar(8, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(6, (1 << 1) + (1 << 0)).writeln();
    gen_squar(6, (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(4, (1 << 1) + (1 << 0)).writeln();
    gen_squar(4, (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(3, (1 << 1) + (1 << 0)).writeln();
    gen_squar(3, (1 << 1) + (1 << 0)).writeln();
    gen_prod_internal(2, (1 << 1) + (1 << 0)).writeln();
    gen_squar(2, (1 << 1) + (1 << 0)).writeln();
}


/** templated finite field F[2^n]
usage:
-----
struct F32
{
    mixin FF!(F32, uint, 32, (1 << 22) + (1 << 2) + (1 << 1) + (1 << 0));
}
-----
Params:
    F = the struct
    uinteger = the storage type
    n = the exponent of two
    p = x^n
*/
mixin template BF(F, uinteger, uint n, ulong p)
{
    mixin declare_prod!(uinteger, n, p);
    mixin (gen_squar(n, p));
    uinteger v;
    this (uinteger v) /// trivial constructor (necessary for implicit conversion from uinteger.)
    {
        this.v = v;
    }
    bool opCast(T)()
    {
        return cast(T)(this.v);
    }

    version (binary)
    string toString() /// toString: binary version
    {
        string ret;
        foreach_reverse (i; 0..n)
        {
            ret ~= ((this.v >> i) & 1).to!string;
        }
        return ret;
    }
    else 
    string toString() /// toString: decimal version
    {
        return this.v.to!string;
    }

    /// equality
    bool opEquals(F other)
    {
        return this.v == other.v;
    }

    /// binary operators with same type
    F opBinary(string op)(F other)
    {
        static if (op == "+" || op == "-")
        {
            return F(this.v ^ other.v);
        }
        else static if (op == "*")
        {
            return F(prod(this.v, other.v));
        }
        else static if (op == "/")
        {
            return this * other.inverse();
        }
    }

    /// ditto
    F opOpAssign(string op)(F other)
    {
        static if (op == "+" || op == "-")
        {
            this.v ^= other.v;
            return this;
        }
        static if (op == "*")
        {
            this.v = prod(this.v, other.v);
            return this;
        }
    }

    /// inverse
    F inverse()
    {
        auto v = squar(this.v);
        auto sq = v;
        foreach (i; 1..(n - 1))
        {
            sq = squar(sq);
            v = prod(v, sq);
        }
        return F(v);
    }

    static if (n == 64)
    {
        private const mask = ulong.max;
    }
    else
    {
        private const mask = (1UL << n) - 1;
    }

    /// power operator
    F opBinary(string op)(ulong e) 
    {
        static if (op == "^^")
        {
            static if (n != 64)
            {
                while (e >> n)
                {
                    e = (e >> n) + (e & mask);
                }
            }
            uinteger v = 1;
            uinteger sq = this.v;
            while (e)
            {
                if (e & 1)
                {
                    v = prod(v, sq);
                }
                sq = squar(sq);
                e >>= 1;
            }
            return F(v);
        }
    }
    F square()
    {
        return F(squar(this.v));
    }
    alias square frobenius;

    bool trace() /// field trace: definition by frobenius automorphism.
    {
        debug (trace) "%s.fieldtrace = ".writef(this);
        F ret;
        F tmp = this;
        foreach (i; 0..n)
        {
            ret += tmp;
            tmp = tmp.square();
        }
        debug (trace) ret.v.writeln();
        return ret.v != 0;
    }
    debug
    {
        private bool another_trace() /// field trace: definition by linear algebra.
        {
            debug (trace) "%s.lineartrace = ".writef(this);
            bool tr;
            foreach (ulong i; 0..n)
            {
                auto before = F(cast(uinteger)(1UL << i));
                if ((this * before).v & before.v)
                {
                    tr = !tr;
                }
            }
            debug (trace) (tr ? 1 : 0).writeln();
            return tr;
        }
        unittest
        {
            version (demo) "checking two definition of trace for x^i: 0 < x < %d ...".writef(n);
            foreach (y; F.standard_basis())
            {
                assert (y.trace == y.another_trace);
            }
            version (demo) " OK".writeln();
        }
    }

    /** field norm over B
    a.norm
        = (a^2^i : 0 <= i < n).prod
        = a ^ (2^i : 0 <= i < n).sum
        = a ^ (2^n-1)
        = 1 (if a != 0); 0 (otherwise)
    */
    bool norm()
    {
        return this.v != 0;
    }

    static F[] standard_basis()
    {
        auto y = F(1), x = F(2);
        F[] ret;
        foreach (i; 0..n)
        {
            ret ~= y;
            y *= x;
        }
        return ret;
    }

    static F trace_one() /// smallest (as normal integer) element which is of trace 1.
    {
        foreach (y; F.standard_basis())
        {
            if (y.trace)
            {
                return y;
            }
        }
        assert (false);
    }
    version (demo) unittest
    {
        "first element with trace 1 in F_{2^%d} is %d".writefln(n, F.trace_one());
    }
    static if (n == 128 || n == 64 || n == 32 || n == 16 || n == 8 || n == 4 || n == 2 || n == 1)
    {
        static F[n] cantor_basis() /// Cantor basis
        out (ret) /// checks linear independency
        {
            uinteger[n] result;
            foreach (i, c; ret)
            {
                result[i] = c.v;
            }
            assert (linear_independency(result));
        }
        body
        {
            auto theta = trace_one();
            F[n] ret;
            ret[0] = F(1);
            foreach (i; 1..n)
            {
                ret[i] = ret[i - 1].next_elem(theta);
            }
            return ret;
        }
        version (demo) unittest
        {
            "Cantor basis: %s".writefln(F.cantor_basis());
        }

        private F next_elem(F theta)
        out (ret)
        {
            auto r = F(ret.v);
            assert (this == r * r + r);
        }
        body
        {
            F ret;
            // [0 <= j < i < n] this^2^j theta^2^i
            uinteger one = 1;
            F thetasq = theta;
            foreach (uinteger i; 0..n)
            {
                F thissq = this;
                foreach (uinteger j; 0..i)
                {
                    ret += thissq * thetasq;
                    thissq = thissq.frobenius();
                }
                thetasq = thetasq.frobenius();
            }
            return ret;
        }
        //static F[n] basis = cantor_basis();
    }
    unittest
    {
        F random_element()
        {
            uinteger v;
            while (v == 0)
            {
                v = random_bits!(uinteger, n)();
            }
            return F(v);
        }
        version (demo) "testing F_{2^%d} ...".writef(n);
        foreach (t; 0..1000)
        {
            F i = random_element(), j = random_element(), k = random_element();
            assert ((i / j) * j == i);
            assert (i * (j+k) == i*j + i*k);
            assert ((i+j) * k == i*k + j*k);
            assert (i * (j*k) == (i*j) * k);
            assert (i / (j*k) == (i/j) / k);
            assert ((i+j) + k == i + (j+k));
        }
        version (demo) " OK".writeln();
    }
}

T random_bits(T, uint n)()
{
    T ret;
    foreach (i; 0..n)
    {
        ret <<= 1;
        if (uniform(0, 2))
        {
            ret |= 1;
        }
    }
    return ret;
}

bool linear_independency(U, size_t n)(U[n] result)
{
outer:
    foreach (i; 0..n)
    {
        foreach (j; 0..n)
        {
            if (result[j] >> i & 1)
            {
                foreach (k; (j+1)..n)
                {
                    if (result[k] >> i & 1)
                    {
                        result[k] ^= result[j];
                    }
                }
                continue outer;
            }
        }
        return false;
    }
    return true;
}
