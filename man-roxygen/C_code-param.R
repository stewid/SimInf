##' @param C_code Character vector with optional model C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded. The DLL
##'     is unloaded and the temporary files are removed after running
##'     the model.
