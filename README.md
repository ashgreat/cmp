# cmp: Conditional Mixed Process Models

[![R-CMD-check](https://github.com/ashgreat/cmp/workflows/R-CMD-check/badge.svg)](https://github.com/ashgreat/cmp/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/cmp)](https://CRAN.R-project.org/package=cmp)

The `cmp` package provides R implementation of Conditional Mixed Process (CMP) models using maximum likelihood estimation. CMP models allow for systems of equations with different model types including continuous, binary and ordered probit, multinomial probit, and censored regression models.

This package is based on David Roodman's Stata CMP package and aims to provide equivalent functionality in R.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ashgreat/cmp")
```

## Quick Start

```r
library(cmp)

# Load example data
data(cmp_example)

# Simple bivariate probit model
model <- cmp(
  formula = list(
    y_continuous ~ x1 + x2 + x3,
    y_binary ~ x1 + x2
  ),
  data = cmp_example,
  indicators = c("continuous", "probit")
)

# View results
summary(model)

# Make predictions
predictions <- predict(model, type = "response")
```

## Supported Model Types

The package supports the following model types through the `indicators` parameter:

- `"continuous"`: Normal linear regression
- `"probit"`: Binary probit model
- `"oprobit"`: Ordered probit model
- `"left"`: Left-censored (Tobit) model
- `"right"`: Right-censored model
- `"interval"`: Interval regression

## Key Features

- **Joint estimation**: Simultaneous estimation of multiple equations with different model types
- **Error correlation**: Accounts for correlation between equation errors
- **Flexible interface**: Easy specification of complex model systems
- **Standard methods**: Comprehensive predict, summary, and diagnostic methods
- **CRAN ready**: Professional package structure with full documentation and testing

## Examples

### Single Equation Models

```r
# Probit model
probit_model <- cmp(
  formula = list(y_binary ~ x1 + x2),
  data = cmp_example,
  indicators = "probit"
)

# Tobit model
tobit_model <- cmp(
  formula = list(y_censored ~ x1 + x2 + x3),
  data = cmp_example,
  indicators = "left"
)
```

### Multiple Equation Systems

```r
# Mixed continuous and binary system
mixed_model <- cmp(
  formula = list(
    y_continuous ~ x1 + x2 + x3,
    y_binary ~ x1 + x2,
    y_ordered ~ x2 + x3
  ),
  data = cmp_example,
  indicators = c("continuous", "probit", "oprobit")
)

summary(mixed_model)
```

## Comparison with Stata CMP

This R implementation provides:

- Compatible model specifications and results
- Maximum likelihood estimation using BFGS optimization
- Support for multiple equation types in single systems
- Proper handling of cross-equation correlation
- Standard model diagnostic and prediction tools

## Documentation

- [Introduction vignette](vignettes/cmp-introduction.Rmd)
- **Function reference:** After installing, run `help(package = "cmp")` for a full index of exported functions.
- **Examples & tutorials:** Review the vignette `vignettes/cmp-introduction.Rmd` or run `vignette("cmp-introduction", package = "cmp")`.

## Citation

If you use this package in your research, please cite:

```
Malshe, Ashwin (2025). cmp: Conditional Mixed Process Models.
R package version 1.0.0. https://github.com/ashgreat/cmp
```

Please also cite the original Stata implementation:

```
Roodman, D. (2011). Fitting fully observed recursive mixed-process models with cmp.
The Stata Journal, 11(2), 159-206.
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

GPL (>= 3)

## Support

- [GitHub Issues](https://github.com/ashgreat/cmp/issues)
- [Package Documentation](https://ashgreat.github.io/cmp/)
