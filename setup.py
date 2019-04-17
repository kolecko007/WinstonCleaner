import setuptools
import os

with open("README.md", "r") as fh:
        long_description = fh.read()

setuptools.setup(
    name="winston_cleaner",
    version="0.1.0",
    author="Serafim Nenarokov, Martin Kolisko",
    author_email="serafim.nenarokov@gmail.com",
    description="WinstonCleaner - transcriptomic data cross-contamination eliminator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kolecko007/WinstonCleaner",
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 2",
                 "Programming Language :: Python :: 2.7",
                 "Operating System :: OS Independent",
    ],
    python_requires='==2.7.*',
    entry_points={
        'console_scripts': ['winston_find_contaminations = bin:find_contaminations.main',
                            'winston_generate_config = bin:generate_config.main',
                            'winston_prepare_data = bin:prepare_data.main'],
    },
    py_modules=['winston'],
    include_package_data=True,
    install_requires=["matplotlib==2.0.2",
                      "numpy>=1.14.0",
                      "scipy>=0.19.0",
                      "biopython>=1.68",
                      "tqdm>=4.14.0",
                      "pathlib>=1.0.1",
                      "PyYAML>=3.11"]
)
