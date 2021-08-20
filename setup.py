import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="go_bench", 
    version="0.0.1",
    author="Andrew Dickson",
    author_email="amdickson@berkeley.edu",
    description="Package for running GO pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/amdson/GO_pipeline",
    packages=["go_bench"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)