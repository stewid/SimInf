## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Create a DESCRIPTION file for the package skeleton
##' @noRd
create_DESCRIPTION_file <- function(path, name, author, maintainer,
                                    email, license)
{
    lines <- c(paste0("Package: ", name),
               "Type: Package",
               paste0("Title: Model ('", name, "') interface to the 'SimInf' package"),
               "Version: 1.0",
               paste0("Author: ", author),
               paste0("Maintainer: ", maintainer, " <", email, ">"),
               paste0("Description: Provides an interface to the 'SimInf' package for the '",
                      name, "' model. More about what it does (maybe more than one line)"),
               paste0("License: ", license),
               paste0("NeedsCompilation: yes"),
               paste0("Depends: SimInf"),
               paste0("LinkingTo: SimInf"))
    writeLines(lines, con = file.path(path, "DESCRIPTION"))
    invisible(NULL)
}
