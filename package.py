# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install rustbca
#
# You can edit this file again by typing:
#
#     spack edit rustbca
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class Rustbca(Package):
    """RustBCA: A Binary Collision Approximation code and libraries for simulating ion-material interactions"""

    homepage = "https://www.github.com/lcpp-org/RustBCA/wiki"
    url      = "https://github.com/lcpp-org/RustBCA/archive/refs/tags/v1.2.0.tar.gz"
    git      = "https://www.github.com/lcpp-org/RustBCA.git"

    # maintainers = ['drobnyjt']

    version('dev', branch='dev')
    version('main', branch='main')
    depends_on('rust')

    def install(self, spec, prefix):
        cargo = which('cargo')
        cargo('build', '--release', '--lib')

        mkdirp(prefix.include)
        install('RustBCA.h', prefix.include)

        mkdirp(prefix.lib)
        install('target/release/liblibRustBCA.so', prefix.lib)
        
