from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="RustBCA",
    version="1.1.2",
    rust_extensions=[
        RustExtension(
            "libRustBCA.pybca",
            binding=Binding.PyO3,
            features=["python"],
            #args=["+nightly", "--edition 2018", "-Z unstable-options"],
            #optional=True,
            rust_version="1.56.1"
            
        )
    ],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)
