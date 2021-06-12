check_program <- function(code) {
    filename <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)

    ## Write the C code to a temporary file.
    c_file <- paste0(filename, ".c")
    on.exit(unlink(c_file), add = TRUE)
    writeLines(code, con = c_file)

    ## The name for the built library from compiling the above C code.
    lib <- paste0(filename, .Platform$dynlib.ext)
    on.exit(unlink(lib), add = TRUE)

    cmd <- c(
        "SHLIB",
        "--clean",
        paste0("--output=", shQuote(lib)),
        shQuote(c_file),
        "-L$(LIB_GSL)/lib",
        "-lm",
        "-lgsl",
        "-lgslcblas")

    invisible(tools::Rcmd(cmd, stderr = NULL, stdout = NULL))
    ifelse(file.exists(lib), "yes", "no")
}
