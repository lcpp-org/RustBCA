from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="RustBCA",
    rust_extensions=[RustExtension("libRustBCA.pybca", binding=Binding.PyO3)],
    # rust extensions are not zip safe, just like C-extensions.
    features=["python"]
    zip_safe=False,
)
