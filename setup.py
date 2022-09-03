from setuptools import setup, find_packages

setup(
    name = "Beaconet",
    version = "1.0.2",
    keywords = ("scRNA-seq","single cell","batch effect","multiple sources","reference-free"),
    description = "A reference-free method for integration of multiple batches of scRNA-seq data.",
    long_description = "the misc-tools for experiment of machine learning and bioinformatics.",
    url = "https://github.com/xuxiaohan",
    author = "Xu Han, Gao Lin",
    author_email = "hxu10670@gmail.com",
    maintainer = "Xu Han",
    maintainer_email = "hxu10670@gmail.com",
    packages = find_packages(),
    include_package_data = True,
    platforms = "any",
    license="MIT Licence",
    install_requires = [
        "numpy",
        "pandas",
        "seaborn",
        "matplotlib",
        "tqdm",
        "scikit-learn",
        "torch",
        "umap-learn",
        "scipy"
        ]
)