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
The subroutine samples a rank 1 array of size \( m \) (`m`) with prescribed 
\( n \) (`n`), the number of elements in the sample. 
It uses the forward Fisher-Yates algorithm to re-shuffle the array (incomplete reshuffle), then generates an index mask (`mask`) from it.

The subroutine only generates the mask rather than a new array of sampled elements. 
Outside the subroutine, the `pack` intrinsic function can simply be used to generate a new array of samples as follows:

`new_array = pack (old_array, mask)`

@note
The procesure has no pure equivalent, because it uses the intrinsic subroutine `random_number` to generate pseudorandom numbers.
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

The subroutine samples a rank 1 array of size \( m \) (`m`). Each element in the array is independently subjected to Bernoulli experiments and included in the sample with specified probability \( m \) (`p`). An index mask (`mask`) for included elements is generated based on these probability experiments.

The subroutine only generates the mask rather than a new array of sampled elements. 
Outside the subroutine, the `pack` intrinsic function can simply be used to generate a new array of samples as follows:

`new_array = pack (old_array, mask)`

@note
The procesure has no pure equivalent, because it uses the intrinsic subroutine `random_number` to generate pseudorandom numbers.
@endnote

### Syntax
`call ` [[fsml(module):fsml_sample_n(interface)]]`(m, p, mask)`

### Parameters
`m`: A scalar of type `integer`.

`p`: A scalar of type `real`. Its value must be between *0.0* and *1.0*.

Invalid argument values will result in returning all `.false.` mask values.

### Returns
`mask`: A rank-1 array of type `logical` with dimension `m`.




<br>
# Examples

```fortran
{!example/example_dat.f90!}
```
