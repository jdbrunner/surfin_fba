import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


try:
    import cplex as cp
except:
    None

try:
    import gurobipy as gb
except:
    None

setuptools.setup(
    name="surfinFBA",
    version="0.8.6",
    author="James D. Brunner, Ph.D.",
    author_email="brunner.james@mayo.edu",
    description="Dynamic FBA for use with COBRApy metabolic models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jdbrunner/surfin_fba",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requres=['numpy','scipy.integrate','pandas','joblib','cobra','matplotlib'],
    include_package_data=True
)
