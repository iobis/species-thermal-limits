# Needed packages on CRAN
req_packages <- c(
    "rstan",
    "dplyr",
    "ggplot2",
    "truncnorm"
)

# Needed packages on GitHub
git_packages <- c("rethinking")
git_packages_source <- c(
    "rmcelreath/rethinking"
)

# Create a function to check if is installed
is_package_installed <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

# Check which ones are not installed and install if needed:
for (i in seq_along(req_packages)) {
    if (!is_package_installed(req_packages[i])) {
        install.packages(req_packages[i])
    }
}

# Check github packages
for (i in seq_along(git_packages)) {
    if (!is_package_installed(git_packages[i])) {
        devtools::install_github(git_packages_source[i])
    }
}
