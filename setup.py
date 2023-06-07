from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="RustBCA",
    version="2.1.0",
    rust_extensions=[
        RustExtension(
            "libRustBCA",
            binding=Binding.PyO3,
            features=["python", "parry3d"],
            #args=["+nightly", "--edition 2021", "-Z unstable-options"],
            #optional=True,
            #rust_version="1.57.0"
        )
    ],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)
