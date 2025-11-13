# Create a package skeleton from a `SimInf_model`

Describe your model in a logical way in R, then `mparse` creates a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object with your model definition that can be installed as an add-on R
package.

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

  The `model`
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object with your model to create the package skeleton from.

- name:

  Character string with the package name. It should contain only (ASCII)
  letters, numbers and dot, have at least two characters and start with
  a letter and not end in a dot. The package name is also used for the
  class name of the model and the directory name of the package.

- path:

  Path to put the package directory in. Default is '.' i.e. the current
  directory.

- author:

  Author of the package.

- email:

  Email of the package maintainer.

- maintainer:

  Maintainer of the package.

- license:

  License of the package. Default is 'GPL-3'.

## Value

invisible `NULL`.

## References

Read the *Writing R Extensions* manual for more details.

Once you have created a *source* package you need to install it: see the
*R Installation and Administration* manual,
[`INSTALL`](https://rdrr.io/r/utils/INSTALL.html) and
[`install.packages`](https://rdrr.io/r/utils/install.packages.html).
