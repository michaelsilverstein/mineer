import setuptools
from mineer import __version__

with open('README.md') as fh:
    long_description = fh.read()

setuptools.setup(
    name='mineer',
    version=__version__,
    entry_points={'console_scripts':[
        'mineer=mineer.pipeline:run',
        'mineer-test-files=mineer.download_test_files:download'
    ]},
    install_requires= ['numpy', 'pandas', 'biopython', 'seaborn', 'numba'],
    url='',
    license='MIT',
    author='michaelsilverstein',
    author_email='michael.silverstein4@gmail.com',
    description='Automated sequence truncation algorithm',
    packages=['mineer'],
    long_description=long_description,
    long_description_content_type = "text/markdown",
    python_requires='>=3.8'
)
