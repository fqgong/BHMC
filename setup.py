import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BHMC-cluster",
    version="0.1.0",
    author="Fuqiang Gong",
    author_email="fqgong@stu.xmu.edu.com",
    description="Small Package to Construct the Structure of Pure Cluster",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fqgong/BHMC",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)"
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy >= 1.19.5",
        "ase >= 3.20.1",
        "dpdata",
        "dpdispatcher",
        "pyxtal"
    ],
    entry_points={
       'console_scripts': [
           'bhmc=BHMC.main:main'
       ]
    }
)
