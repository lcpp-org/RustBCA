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
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    homepage = "https://www.github.com/lcpp-org/RustBCA/wiki"
    url      = "https://github.com/lcpp-org/RustBCA/archive/refs/tags/v1.0.0.tar.gz"
    git      = "https://www.github.com/lcpp-org/RustBCA.git"

    # FIXME: Add a list of GitHub accounts to
    # notify when the package is updated.
    # maintainers = ['github_user1', 'github_user2']

    version('dev', branch='dev')
    version('1.0.0', sha256='99dcac7c7a78e6cd17da63a0dcbb3c36bca523ffafbb0425128b0c971b1a6829')
    depends_on('rust')

    # FIXME: Add dependencies if required.
    # depends_on('foo')

    def install(self, spec, prefix):
        cargo = which('cargo')
        cargo('build', '--release', '--lib', '--target-dir', prefix)
