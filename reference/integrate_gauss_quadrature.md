# Gaussian quadrature for mixture density integration

Provides S3 methods for `integrate_density` and `integrate_density_log`
that use Gaussian quadrature matched to the mixture family:

- `normMix`: Gauss-Hermite quadrature

- `betaMix`: Gauss-Jacobi quadrature

- `gammaMix`: Gauss-Laguerre quadrature

## Details

The integration method is controlled by the global option
`RBesT.integrate_method` (default `"GQ"`). Setting it to `"adaptive"`
falls through to the default adaptive method.

Tolerance control is enabled by default: starting from `RBesT.GQ_nodes`
nodes, the node count is repeatedly increased (by
`RBesT.GQ_node_growth`, default doubling) until successive estimates
agree within `max(RBesT.GQ_abs_tol, RBesT.GQ_rel_tol * |I|)`. Setting
`RBesT.GQ_rel_tol` to a non-finite or non-positive value (e.g. `Inf`)
restores the legacy single-shot evaluation at `RBesT.GQ_nodes` nodes.

Relevant options (all consulted only on the `"GQ"` path):

- `RBesT.GQ_nodes` (default `20L`) – starting node count

- `RBesT.GQ_rel_tol` (default `1e-4`) – relative tolerance; a non-finite
  or non-positive value (e.g. `Inf`) disables refinement (single
  evaluation at `RBesT.GQ_nodes`)

- `RBesT.GQ_abs_tol` (default `1e-6`) – absolute tolerance floor

- `RBesT.GQ_max_nodes` (default `240L`) – refinement cap

- `RBesT.GQ_node_growth` (default `2`) – node growth factor

- `RBesT.GQ_on_nonconvergence` (default `"adaptive"`) – one of
  `"adaptive"` (fall through to adaptive integration), `"warn"`,
  `"error"`, `"silent"`
