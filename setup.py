from setuptools import setup
from mineer import __version__

setup(
    name='mineer',
    version=__version__,
    entry_points={'console_scripts':['mineer=mineer.pipeline:mineer_cli']},
    install_requires= ['numpy', 'pandas', 'biopython'],
    url='',
    license='MIT',
    author='michaelsilverstein',
    author_email='michael.silverstein4@gmail.com',
    description='Automated sequence truncation algorithm',
    long_description_content_type = "text/markdown",
)
