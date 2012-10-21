module fftu;

import std.stdio;

import metafft;
import fint;

mixin template makeFFT(F, U, uint n)
{
    mixin FFTTOOLS!(F, U) mytools;
    struct FFT
    {
        mixin WZC!(F, U, n);
    }
    FFT fft;
    unittest
    {
        mytools.fftprod([10], [17]).writeln();
    }
}

mixin makeFFT!(B4, ubyte, 4) fft4;
mixin makeFFT!(B8, ubyte, 8) fft8;
mixin makeFFT!(B16, ushort, 16) fft16;
mixin makeFFT!(B32, uint, 32) fft32;
