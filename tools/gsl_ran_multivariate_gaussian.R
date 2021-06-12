source("tools/check_program.R")

code <- c(
    "#include <gsl/gsl_randist.h>",
    "int main()",
    "{",
    "    gsl_ran_multivariate_gaussian(NULL, NULL, NULL, NULL);",
    "    return 0;",
    "}")

cat(check_program(code))
