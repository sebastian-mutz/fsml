---
title: DAT: Data Manipulation
---

# Overview

This is the API documentation for data manipulation and processing procedures.

[TOC]

<br>
# Fixed-N Sampling

## `fsml_sample_n`

### Description
The subroutine samples a rank-1 array of size \( m \) (`m`) with prescribed 
\( n \) (`n`), the number of elements in the sample. It uses the forward 
Fisher-Yates algorithm to re-shuffle the array (incomplete reshuffle), 
then generates an index mask (`mask`) from it.

The subroutine only generates the mask rather than a new array of sampled
elements. Outside the subroutine, the `pack` intrinsic function can simply 
be used to generate a new array of samples as follows:

```Fortran
new_array = pack(old_array, mask)
```

@note
The procesure has no pure equivalent, because it uses the intrinsic
subroutine `random_number` to generate pseudorandom numbers.
@endnote


### Syntax
`call ` [[fsml(module):fsml_sample_n(interface)]]`(m, n, mask)`

### Parameters
`m`: A scalar of type `integer`.

`n`: A scalar of type `integer`. It must be non-zero positive and not larger than `m`.

Invalid argument values will result in returning all `.false.` mask values.

### Returns
`mask`: A rank-1 array of type `logical` with dimension `m`.


<br>
# Probabilistic (Poisson) Sampling

## `fsml_sample_p`

### Description

The subroutine samples a rank-1 array of size \( m \) (`m`). Each element
in the array is independently subjected to Bernoulli experiments to determine
its inclusion or exclusion from the sample. The inclusion probability is specified
through \( p \) (`p`). An index mask (`mask`) for included elements is generated
based on these probability experiments.

The subroutine only generates the mask rather than a new array of sampled
elements. Outside the subroutine, the `pack` intrinsic function can simply 
be used to generate a new array of samples as follows:

```Fortran
new_array = pack(old_array, mask)
```

@note
The procesure has no pure equivalent, because it uses the intrinsic
subroutine `random_number` to generate pseudorandom numbers.
@endnote

### Syntax
`call ` [[fsml(module):fsml_sample_p(interface)]]`(m, p, mask)`

### Parameters
`m`: A scalar of type `integer`.

`p`: A scalar of type `real`. Its value must be between *0.0* and *1.0*.

Invalid argument values will result in returning all `.false.` mask values.

### Returns
`mask`: A rank-1 array of type `logical` with dimension `m`.


<br>
# K-Fold Sampling

## `fsml_sample_k`

### Description

The subroutine samples a rank-1 array of size \( m \) (`m`). It creates
\( k \) (`k`) ~equal-sized samples of the rank-1 array. The array
indices are shuffled using the Fisher–Yates algorithm. Then, \( k \) logical
masks are constructed. In each mask, the indices belonging to one of the 
\( k \) folds (not part of the sample) are set to `.false.` and the remaining
indices are set to `.true.`, making the masks directly suitable for
k-fold cross-validation.

In cases where \( m \) cannot be divided into exactly equal sized \( k \) folds,
the remainder is added to the last fold. Therefore, the last sample may be smaller
than the others.

The masks can be applied using the `pack` intrinsic function. For the \( i^{th} \) of the \( k \) samples, this could be done as follows:

```Fortran
new_array = pack(old_array, mask(:, i)) 
```

@note
The procesure has no pure equivalent, because it uses the intrinsic
subroutine `random_number` to generate pseudorandom numbers.
@endnote

### Syntax
`call ` [[fsml(module):fsml_sample_k(interface)]]`(m, k, mask)`

### Parameters
`m`: A scalar of type `integer`.

`k`: A scalar of type `integer`. It must be non-zero positive and not larger than `m`.

Invalid argument values will result in returning all `.false.` mask values.

### Returns
`mask`: A rank-2 array of type `logical` with dimensions `m`, `k`.


<br>
# Examples

```fortran
{!example/example_dat.f90!}
```
