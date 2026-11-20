# RBesT: R Bayesian evidence synthesis tools

Tool-set to support Bayesian evidence synthesis. This includes
meta-analysis, (robust) prior derivation from historical data,
operating characteristics and analysis (1 and 2 sample cases). Please
refer to [Weber et al. (2021)](https://doi.org/10.18637/jss.v100.i19)
for details on applying this package while [Neuenschwander et
al. (2010)](https://doi.org/10.1177/1740774509356002) and [Schmidli et
al. (2014)](https://doi.org/10.1111/biom.12242) explain details on the
methodology.


## Getting Started

To get a quick overview over the package functionality, please refer
to the introductory [vignettes](https://opensource.nibr.com/RBesT/articles/index.html). We recommend
starting with the [Getting started with `RBesT`
(binary)](https://opensource.nibr.com/RBesT/articles/introduction.html) vignette as a first read.

For a brief 15 minute introduction to using `RBesT`, have a look at
the [ESMARConf 2023
tutoral](https://youtu.be/w90j-PLwWDc?si=coN1t1emYPGlT5mT).

The `rbest-assistant` AI skill helps assistants such as Claude, Posit
Assistant, GitHub Copilot or ChatGPT Codex guide RBesT workflows, from
MAP priors to operating characteristics. See the [AI assistant
skill](https://opensource.nibr.com/RBesT/articles/rbest-assistant.html)
article for installation. Then try:

> Using the RBesT assistant skill: give me a short introduction to MAP
> priors. How can I use historical control data for a new trial? No
> analysis yet, just an overview. Then briefly tell me what else you can
> help me with in RBesT.

The agent's output is not validated by `RBesT` and can differ between
runs; check the code and results, and have any analysis or design
reviewed by a qualified statistician. The assistant is an experimental
feature and may change; watch for updates, and we welcome your feedback
via [GitHub issues](https://github.com/Novartis/RBesT/issues).
