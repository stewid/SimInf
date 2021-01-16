VERSION <- commandArgs(TRUE)
if (!file.exists(sprintf("../windows/gsl-%s/include/gsl/gsl_rng.h", VERSION))) {
    url <- sprintf("https://github.com/rwinlib/gsl/archive/v%s.zip", VERSION)
    download.file(url, "lib.zip", quiet = TRUE)
    dir.create("../windows", showWarnings = FALSE)
    unzip("lib.zip", exdir = "../windows")
    unlink("lib.zip")
}
