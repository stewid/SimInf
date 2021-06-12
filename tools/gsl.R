source("tools/check_program.R")

code <- c(
    "#include <gsl/gsl_rng.h>",
    "int main()",
    "{",
    "    gsl_rng *r = NULL;",
    "    return 0;",
    "}")

cat(check_program(code))
