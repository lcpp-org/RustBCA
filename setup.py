from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="RustBCA",
    version="3.0.0",
    rust_extensions=[
        RustExtension(
            "libRustBCA",
            binding=Binding.PyO3,
            features=["python", "parry3d", "pythonize", "cpr_rootfinder"],
        )
    ],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)
