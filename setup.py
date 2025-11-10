# https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/
# https://packaging.python.org/en/latest/tutorials/packaging-projects/

from setuptools import setup, find_packages

VERSION = '0.0.0a1'
DESCRIPTION = 'Tulips3D data prep'
LONG_DESCRIPTION = ''

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="DataPrepTulips3D", 
        version=VERSION,
        author="B.L. de Vries",
        author_email="<bldevries@protonmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=["numpy>=2.0.2", "scipy>=1.13.1"], # add any additional packages that 
        license = " GPL-3.0",
        # needs to be installed along with your package. Eg: 'caer'
        url= "https://github.com/bldevries/tulips3D",
        keywords=['python'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)