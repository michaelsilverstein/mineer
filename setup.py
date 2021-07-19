import setuptools
from minEER import __version__

setuptools.setup(
    name='mineer',
    version=__version__,
    entry_points={'console_scripts':['mineer=mineer.pipeline:mineer_cli']},
    install_requires= ['numpy', 'pandas', 'biopython'],
    url='',
    license='MIT',
    author='michaelsilverstein',
    author_email='michael.silverstein4@gmail.com',
    description='Automated sequence truncation algorithm',
    packages=setuptools.find_packages(where="mineer"),
    long_description_content_type = "text/markdown",
)
