# Contributing to the Fortran Statistics and Machine Learning library (FSML)

Thank you for considering contributing to FSML. Please review and follow the guidelines below to make the process
efficient and enjoyable for all. Since this is a relatively small project, the guidelines are kept short and you
are encouraged to simply contact [sebastian@mutz.science](mailto:sebastian@mutz.science) when in doubt.

If you are interested in contributing code, you certify that you own or are allowed to share the content of your
contribution under the [fsml license](https://github.com/sebastian-mutz/fsml/blob/HEAD/LICENCE).

## Style

A consistent style improves the readability and maintainability of software. The best way to get a sense for it
is to study the existing code. Some of the conventions for the project are highlighted below:

### Use modern standard Fortran and keep it simple!

* We have thankfully moved on from `goto` statements and the likes. Best leave them in the past. :)
* Do not use vendor extensions or non-standard syntax like `real*8`.
* Use basic constructs if they solve the problem. Keeping it simple makes the code more readable and maintainable.

### Indentation, spaces, and line length

* The body of every Fortran construct should be indented by **two spaces**.
* A line should not exceed **80 characters**; use line breakers for long expressions.
* Use **spaces**, not tabs for indentation.
* Use **3 spaces of indentation** for `if`, `do`, etc..

### Variable and procedure names

* Variable names should be made up understandable lowercase words or abbreviations linked by underscore
  (e.g., `eof_analysis`).
* Function and subroutine should follow the naming conventions `f_[module_name]_[function_descriptor]`
  and `s_[module_name]_[function_descriptor]`, respectively. Pure procedures that have impure wrappers
  should append `_core`. For example, `f_nlp_kmeans_core`. Public-facing interfaces should simplify this
  to `fsml_kmeans`.
* Append the `program`, `module`, `function` or `subroutine` names in their end statement.
  For example, `end subroutine f_nlp_kmeans_core`.
  
### Kinds and Attributes

* All `real` should use the working precision (`real(wp)`), and integers should be `integer(i4)`, unless
  there is a good reason to deviate from it.
* Always specify the `intent` for arguments. The `optional` attribute should follow the `intent`.
* Don't use the `dimension` attribute. Instead, use  `real, allocatable :: x(:)`

### Documentation

* FORD documentation strings should be provided for all public modules, procedures, arguments, and internal variables.
* Documentation for every new procedure should immediately be added to the modules respective documentation file
  (see files in `doc/api`)


## Opening an Issue (Bug Reports and Suggestions)

### Bug Report

These are extremely valuable for improving software. If you can identify a problem
with the code, follow these simple steps in reporting it:

* Make sure you have the latest version of the software.
* Check if it has already been reported.the issue has already been reported.
* Create a minimal example that reproduces the problem.
* In your report, mention the version and compiler you are using and include the minimal example,
  describe what you expect the code to do and what it actually does.
  
### Feature Suggestion

* Think about whether or not the feature fits the scope and aims of the software and check if it has already been discussed.
* If it is a new idea that fits the scope and aims, explain the feature concisely and make an argument for it.
  Why is it needed? How does it support the software's aims? Examples of similar implementations help to discuss it.

