ClonalIndices
================

# ClonalIndices

The `ClonalIndices` package provides computational tools for clonal
selection based on predicted genotypic effects. It implements a set of
selection indices used to rank genotypes according to multiple trait
information.

## Key Features

The package currently includes the following indices:

- Balanced Selection Index (BSI) — newly developed index

- Multiplicative Balanced Selection Index (MBSI) — newly developed index

- Adapted Smith Selection Index (SSI) — (Smith, 1936)

- Adapted Elston Selection Index (ESI) — (Elston, 1963)

- Adapted Mulamba & Mock Selection Index (MSI) — (Mulamba & Mock, 1978)

## Dataset

The ClonalIndices package includes a dataset, `Arinto`, used to
demonstrate the application of the selection indices implemented in the
package.

The dataset contains predicted genotypic information for a set of clones
and is provided as an example to illustrate the workflow of index
computation and genotype ranking.

It is used in the example of the main functions to ensure full
reproducibility of the methodology presented in Surgy et al (in
preparation).

## References

Elston, R., 1963. A weight-free index for the purpose of ranking or
selection with respect to several traits at a time. Biometrics, Volume
19, pp. 85-97.

Mulamba, N. N. and Mock, J., 1978. Improvement of yield potential of the
ETO blanco maize (Zea mays L.) population by breeding for plant traits.
Egyptian Journal of Genetics and Cytology, Volume 7, pp. 40-57.

Smith, H., 1936. A discriminant function for plant selection. Ann.
Eugen, Volume 7, pp. 240-250.

Surgy, S., Cadima, J., Gonçalves, E.. New Selection Indices for
Selection within ancient grapevine varieties (in preparation).

Surgy, S., Pinto, T., Faísca-Silva, D., Coelho, J., Martins, A., &
Gonçalves, E., 2024. Multi-trait selection in ancient grapevine
varieties. IVES Conference Series.

Surgy, S., 2026. New methods and tools for multi-trait selection in
grapevine breeding programmes. Doctoral thesis submitted to the
Instituto Superior de Agronomia, Universidade de Lisboa (pending
defence).

## Methodological framework

The theoretical framework, including the normalization procedure and
biological interpretation of these indices, is currently under
preparation for publication. The following sections present the
definition of each index as implemented in the package.

### Normalization of Empirical Best Linear Unbiased Predictors (EBLUPs) of genotypic effects

The Empirical Best Linear Unbiased Predictors (EBLUPs ) of genotypic
effects resulting from the fitting of appropriate mixed models were
normalized as follows:

$$
EBLUP\%_{k_i} =
\frac{EBLUP_{k_i}}{\overline{y}_{k_{..}}}
\times 100
$$

where:

- $EBLUP\%_{k_i}$ is the normalized EBLUP of genotypic effect of clone
  *i* for trait *k*  
- $EBLUP_{k_i}$ is the EBLUP of genotypic effect of clone *i* for trait
  *k*  
- $\overline{y}_{k_{..}}$ is the overall mean of observations for trait
  *k*

Normalized EBLUPs of genotypic effects are used as the basis for the
computation of BSI, SSI and MSI selection indices.

### Balanced selection index (BSI)

For a given clone *i*, the BSI is defined as:

$$
BSI_i = \frac{\sum_{k=1}^{m} EBLUP\%_{k_i}}{\sum_{k=1}^{m} \left|EBLUP\%_{k_i}\right|}
$$

where:

- $EBLUP\%_{k_i}$ is the normalized EBLUP of genotypic effect of clone
  *i* for trait *k*  
- $m$ is the total number of traits considered

### Multiplicative Balanced selection index (MBSI)

For a given clone *i*, the MBSI is defined as:

$$
\mathrm{MBSI}_{i}=
\prod_{k=1}^{m}
\left(
\left|EBLUP\%_{k_i}\right|
+1
\right)^{r}
$$ where:

- $\left|EBLUP\%_{k_i}\right|$ is the absolute value of the normalized
  EBLUP of genotypic effect of clone *i* for trait *k*  
- $m$ is the total number of traits considered  
- $r = 1$ when ${EBLUP}_{k_{i}}\geq 0$ and $r = -1$ when
  ${EBLUP}_{k_{i}} \lt 0$.

### Adapted Smith selection index (SSI)

The SSI for clone *i* is defined as:

$$
SSI_i = \sum_{k=1}^{m} EBLUP\%_{k_i}
$$

where:

- $EBLUP\%_{k_i}$ is the normalized EBLUP of genotypic effect of clone
  *i* for trait *k*  
- $m$ is the total number of traits considered

### Adapted Elston selection index (ESI)

The ESI for clone *i* is defined as:

$$
ESI_i = \prod_{k=1}^{m} \left( (p_{k_i} - p_{b_k}) \cdot \delta_k \right)
$$

where:

- $p_{k_i}$ is the predicted genotypic value obtained by the fitting of
  linear mixed models, of clone *i* for trait *k*
- $p_{b_k}$ is the baseline value for trait *k*, defined as the least
  desirable predicted genotypic value for that trait  
- $\delta_k$ is a directional coefficient such that:
  - $\delta_k = -1$ for traits to be decreased  
  - $\delta_k = 1$ for traits to be increased  
- $m$ is the total number of traits considered

### Adapted Mulamba & Mock selection index (MSI)

The MSI for clone *i* is defined as:

$$
MSI_i = \sum_{k=1}^{m} r_{k_i}
$$

where:

- $r_{k_i}$ is the rank of clone *i* for trait *k*, based on normalized
  EBLUP of genotypic effect and adjusted according to the selection
  direction  
- $m$ is the total number of traits considered

## Installation

To install from CRAN:

``` r
install.packages("ClonalIndices")
```

## Dependencies

This package depends on the following R packages:

- **stats**: base R package used for common statistical and utility
  functions.

## Example

In this example, we ask for the top 10 clones in the ranking obtained
with SSI. With the object `rankall` all the indices are computed if
values are indicated for it (EBLUPs of genotypic effects for BSI, MBSI,
SSI and MSI, and predicted genotypic values for ESI).

``` r
library(ClonalIndices)
indexcl <- clonind(
          blups = c("ydblup", "bwblup"),
          preds = c("ydpred", "bwpred"),
          clnum = 10,
          meanvec = c(ydblup = 4.841, bwblup = 1.342),
          criteria = c(ydblup = 1, bwblup = -1, ydpred = 1,bwpred = -1),
          indextype = "ssi",
          data = Arinto
          )

indexcl
#> Rank of genotypes for index SSI 
#>      Ref      SSI
#> 1  AR077 48.94075
#> 2  AR073 35.64610
#> 3  AR069 30.97629
#> 4  AR083 27.66782
#> 5  AR081 27.44441
#> 6  AR112 25.39052
#> 7  AR054 25.35562
#> 8  AR131 24.32555
#> 9  AR080 24.27698
#> 10 AR015 22.94989

indexcl$rankall
#>    Ref.BSI BSI Ref.MBSI      MBSI Ref.SSI      SSI Ref.MSI MSI Ref.ESI
#> 1    AR083   1    AR083 108.36705   AR077 48.94075   AR083 303   AR083
#> 2    AR081   1    AR161 105.34685   AR073 35.64610   AR131 297   AR161
#> 3    AR112   1    AR131  97.19992   AR069 30.97629   AR161 292   AR131
#> 4    AR131   1    AR081  79.74767   AR083 27.66782   AR081 290   AR081
#> 5    AR080   1    AR093  70.32264   AR081 27.44441   AR020 278   AR093
#> 6    AR015   1    AR114  70.29367   AR112 25.39052   AR100 275   AR112
#> 7    AR003   1    AR112  68.63031   AR054 25.35562   AR112 274   AR114
#> 8    AR020   1    AR020  64.20859   AR131 24.32555   AR114 272   AR020
#> 9    AR161   1    AR100  58.93775   AR080 24.27698   AR093 270   AR117
#> 10   AR023   1    AR117  58.34819   AR015 22.94989   AR117 263   AR130
#>          ESI
#> 1  0.8320670
#> 2  0.8186860
#> 3  0.8012759
#> 4  0.7816492
#> 5  0.7523872
#> 6  0.7519009
#> 7  0.7384095
#> 8  0.7294277
#> 9  0.7285551
#> 10 0.7066715
```

## Citation

This package can be cited as follows:

Surgy, S., Cadima, J. & Gonçalves, E., 2026. ClonalIndices - Clonal
Selection Indices Based on Predicted Genotypic Effects.

## License

This package is released under the GPL-3 License.

## Contact and Bug reports

For suggestions, and bug reports use the link:

https://github.com/soniasurgy/ClonalIndices/issues

