OptErd
======

The code can be build with [Ninja build system](http://martine.github.io/ninja/). The script configure.py generate Ninja build file (build.ninja) which builds the libraries (ERD and OED), C interface, and tests.

Requirements
------------

* Python 2 with `ninja` and `argparse` packages
* [Ninja build system](http://martine.github.io/ninja/)
* Intel Compiler
* Intel MPSS for MIC builds (including builds with MIC offload)

Configuration

To generate Ninja build files execute `./configure.py <options...>`. Try `./configure --help` for the list of options.

Build
-----

After you configured the build, go to the root `OptErd` directory and run `ninja` to build the project. The static libraries will be built in `lib` and the tests will go to `test` directory.

If you want to remove generated files, run `ninja -t clean` (**note**: usually it is not needed as Ninja tracks all dependencies between files, and if you change a C header, all source files that include will be recompiled automatically).
