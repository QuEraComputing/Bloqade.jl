# ComplexArrays

Struct of array layout for complex number.

`StructArrays` implements a similar layout on
discontinuous memory, thus it doesn't work with
BLAS well. This package is mainly aimming to
provide a specialized array for complex array
that works with `Adapt` interface.
