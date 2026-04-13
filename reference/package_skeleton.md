# Create a package skeleton from a `SimInf_model`

Generate a source package directory from a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object, allowing the model to be installed as a standalone add-on R
package. This is useful for distributing custom models or for improving
startup performance by compiling the C code once during package
installation rather than on the first call to
[`run()`](http://stewid.github.io/SimInf/reference/run.md).

## Usage

``` r
package_skeleton(
  model,
  name = NULL,
  path = ".",
  author = NULL,
  email = NULL,
  maintainer = NULL,
  license = "GPL-3"
)
```

## Arguments

- model:

  A
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object to create the package skeleton from.

- name:

  Character string with the package name. It must contain only ASCII
  letters, numbers, and dots; have at least two characters; start with a
  letter; and not end in a dot. The package name is also used for the
  class name of the model and the directory name of the package.

- path:

  Path to put the package directory in. Default is `"."` (the current
  directory).

- author:

  Author of the package.

- email:

  Email address of the package maintainer.

- maintainer:

  Maintainer of the package.

- license:

  License of the package. Default is `"GPL-3"`.

## Value

`invisible(NULL)`.

## Details

The typical workflow is:

1.  Define the model using
    [`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) or
    [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md).

2.  Call `package_skeleton` to generate the package

3.  Install the package using
    [`install.packages`](https://rdrr.io/r/utils/install.packages.html)
    or `R CMD INSTALL`.

## References

Read the *Writing R Extensions* manual for more details on building R
packages.

## See also

[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) for
creating a `SimInf_model` object from a model specification.
[`INSTALL`](https://rdrr.io/r/utils/INSTALL.html) and
[`install.packages`](https://rdrr.io/r/utils/install.packages.html) for
installing the generated package.
