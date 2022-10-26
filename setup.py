from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.1'
DESCRIPTION = 'A Package that renders different Calculations in Latex based on the SIA-Norms used by Structural Engineers.'
LONG_DESCRIPTION = 'Currently there are plenty of Excel-sheets used to do Calculations in the structural Engineering-field. This Package includes the most used Verification-Calculations....tbc'

# Setting up
setup(
    name="Handnachweise",
    version=VERSION,
    author="Pascal Gitz",
    author_email="<pascal.gitz@hotmail.ch>",
    description=DESCRIPTION,
    packages=find_packages(),
    install_requires=[],
    keywords=['handcalcs', 'structural Engineer', ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Structural Engineers from Switzerland",
        "Programming Language :: Python :: 3",
        "Operating System :: Microsoft :: Windows",
    ]
)

