import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="surfin_fba", # Replace with your own username
    version="0.0.1",
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
)
