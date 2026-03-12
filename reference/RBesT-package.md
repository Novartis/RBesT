# R Bayesian Evidence Synthesis Tools

The RBesT tools are designed to support in the derivation of parametric
informative priors, asses design characeristics and perform analyses.
Supported endpoints include normal, binary and Poisson.

## Details

For introductory material, please refer to the vignettes which include

- Introduction (binary)

- Introduction (normal)

- Customizing RBesT Plots

- Robust MAP, advanced usage

The main function of the package is
[`gMAP()`](https://opensource.nibr.com/RBesT/reference/gMAP.md). See
it's help page for a detailed description of the statistical model.

## Global Options

|                            |                                     |                                                                                     |
|----------------------------|-------------------------------------|-------------------------------------------------------------------------------------|
| Option                     | Default                             | Description                                                                         |
| `RBesT.MC.warmup`          | 2000                                | MCMC warmup iterations                                                              |
| `RBesT.MC.iter`            | 6000                                | total MCMC iterations                                                               |
| `RBesT.MC.chains`          | 4                                   | MCMC chains                                                                         |
| `RBesT.MC.thin`            | 4                                   | MCMC thinning                                                                       |
| `RBesT.MC.save_warmup`     | `FALSE`                             | MCMC warmup samples saving                                                          |
| `RBesT.MC.control`         | `list(adapt_delta=0.99,`            | sets `control` argument for Stan call                                               |
|                            | `stepsize=0.01,`                    |                                                                                     |
|                            | `max_treedepth=20)`                 |                                                                                     |
| `RBesT.MC.ncp`             | 1                                   | parametrization: 0=CP, 1=NCP, 2=Automatic                                           |
| `RBesT.MC.init`            | 1                                   | range of initial uniform \\\[-1,1\]\\ is the default                                |
| `RBesT.MC.rescale`         | `TRUE`                              | Automatic rescaling of raw parameters                                               |
| `RBesT.verbose`            | `FALSE`                             | requests outputs to be more verbose                                                 |
| `RBesT.integrate_args`     | `list(lower=-Inf,`                  | arguments passed to `integrate` for                                                 |
|                            | `upper=Inf,`                        | intergation of densities                                                            |
|                            | `rel.tol=.Machine$double.eps^0.25,` |                                                                                     |
|                            | `abs.tol=.Machine$double.eps^0.25,` |                                                                                     |
|                            | `subdivisions=1E3)`                 |                                                                                     |
| `RBesT.integrate_prob_eps` | `1E-6`                              | probability mass left out from tails if integration needs to be restricted in range |

## Version History

See `NEWS.md` file.

## References

Stan Development Team (2020). RStan: the R interface to Stan. R package
version 2.19.3. https://mc-stan.org

## See also

Useful links:

- <https://opensource.nibr.com/RBesT/>

- Report bugs at <https://github.com/Novartis/RBesT/issues>

## Author

**Maintainer**: Sebastian Weber <sebastian.weber@novartis.com>

Other contributors:

- Novartis Pharma AG \[copyright holder\]

- Beat Neuenschwander <beat.neuenschwander@novartis.com> \[contributor\]

- Heinz Schmidli <heinz.schmidli@novartis.com> \[contributor\]

- Baldur Magnusson <baldur.magnusson@novartis.com> \[contributor\]

- Yue Li <yue-1.li@novartis.com> \[contributor\]

- Satrajit Roychoudhury <satrajit.roychoudhury@novartis.com>
  \[contributor\]

- Lukas A. Widmer <lukas_andreas.widmer@novartis.com>
  ([ORCID](https://orcid.org/0000-0003-1471-3493)) \[contributor\]

- Daniel Sabanés Bové <daniel.sabanes_bove@rconis.com>
  ([ORCID](https://orcid.org/0000-0002-0176-9239)) \[contributor\]

- Trustees of Columbia University (R/stanmodels.R, configure,
  configure.win) \[copyright holder\]
