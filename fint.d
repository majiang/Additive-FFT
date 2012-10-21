module fint;

import std.stdio;
import std.conv : to;


import metafint;

version = large;
version = small;

version (small)
{
    struct B2
    {
        mixin BF!(B2, ubyte, 2, (1 << 1) + (1 << 0));
    }
    struct B3
    {
        mixin BF!(B3, ubyte, 3, (1 << 1) + (1 << 0));
    }
    struct B4
    {
        mixin BF!(B4, ubyte, 4, (1 << 1) + (1 << 0));
    }
    struct B6
    {
        mixin BF!(B6, ubyte, 6, (1 << 1) + (1 << 0));
    }
    struct B8
    {
        mixin BF!(B8, ubyte, 8, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0));
    }
    struct B12
    {
        mixin BF!(B12, ushort, 12, (1 << 10) + (1 << 2) + (1 << 1) + (1 << 0));
    }
    struct B16
    {
        mixin BF!(B16, ushort, 16, (1 << 12) + (1 << 3) + (1 << 1) + (1 << 0));
    }
}
struct B24
{
    mixin BF!(B24, uint, 24, (1 << 7) + (1 << 2) + (1 << 1) + (1 << 0));
}
struct B32
{
    mixin BF!(B32, uint, 32, (1 << 22) + (1 << 2) + (1 << 1) + (1 << 0));
}
version(large) struct B64
{
    mixin BF!(B64, ulong, 64, (1UL << 61) + (1UL << 34) + (1UL << 9) + (1UL << 0));
}
