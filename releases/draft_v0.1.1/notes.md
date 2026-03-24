# FSML v0.1.1
 
## Changelog
The following has changed since the last release (v0.1.0):

### Bug Fixes
- Fixed gamma distribution loc bug and removed redundant arg checks in pure procedure ([f578d27](https://github.com/sebastian-mutz/fsml/commit/f578d27fba2518bb9b0f66056f2a59fb233fbf25)) [contributed by @sebastian-mutz]

### New Features

#### LIN module
- [ ] Added Least Absolute Shrinkage and Selection Operator (LASSO) regression [contributed by @loiseaujc].
- [ ] Added Canonical Correlation Analysis (CCA) [contributed by @loiseaujc].

#### STS module
- [x] Added function to compute quantile/percentile value (based on Hyndman & Fan (1996) type 7 index). ([562d9af](https://github.com/sebastian-mutz/fsml/commit/562d9af29dccff7c25de1766039ecd304ffc33c4), [a0c67ef](https://github.com/sebastian-mutz/fsml/commit/a0c67ef4b85275220860381db073199b33797884), [262f21b](https://github.com/sebastian-mutz/fsml/commit/262f21b31b47b06267d06a2e43899eed0cc4858b)) [contributed by @sebastian-mutz]

#### DST module
- [x] Added logistic distribution PDF, CDF, and PPF. ([3f81edb](https://github.com/sebastian-mutz/fsml/commit/3f81edb9e9ecfe48a34a2d4a75115fb909d9cccc), [3ae5c64](https://github.com/sebastian-mutz/fsml/commit/3ae5c6413967e13ff11495414b4d83f26f103453), [a287f34](https://github.com/sebastian-mutz/fsml/commit/a287f34ca90b05c863bb840351e641e5f36cab02), [080bf3b](https://github.com/sebastian-mutz/fsml/commit/080bf3b2f8cc6da3158bcdc05418d3c9792835c7)) [contributed by @sebastian-mutz]
- [x] Added log-logistic distribution PDF, CDF, and PPF. ([3ad09b9](https://github.com/sebastian-mutz/fsml/commit/3ad09b9107c1aaa6f5d0516491d32a1ac3d7999), [3ae5c64](https://github.com/sebastian-mutz/fsml/commit/3ae5c6413967e13ff11495414b4d83f26f103453), [a287f34](https://github.com/sebastian-mutz/fsml/commit/a287f34ca90b05c863bb840351e641e5f36cab02), [080bf3b](https://github.com/sebastian-mutz/fsml/commit/080bf3b2f8cc6da3158bcdc05418d3c9792835c7)) [contributed by @sebastian-mutz]

#### NLP module

- [ ] Implemented Lance-Williams algorithm for clustering; more efficient than re-computing Mahalanobis distance [contributed by @sebastian-mutz]

## New Contributors
- @loiseaujc/[Jean-Christophe Loiseau](https://loiseaujc.github.io/) (PRs).

## Communication

- [ ] Conference Presentation: FSML's development approach was presented at the EGU 2026 General Assembly session ["Scientific software development in the Geosciences: good practices, pitfalls and solutions"](https://www.egu26.eu/session/57206) to encourage more sustainable development practices in the scientific community. [contributed by @sebastian-mutz]
- [ ] Blog post: This blog post is a brief development update, providing more context and background. [contributed by @sebastian-mutz]

## Application Spotlight
- @Khadar146/Khadar Dahir Abdisalan used FSML's distribution functions to statistically model future droughts in Somaliland.

