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
    name="surfin_fba",
    version="0.6",
    author="James Brunner, Ph.D.",
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
    install_requres=['numpy','scipy.integrate','pandas','time','joblib','cobra']
)
