from setuptools import setup
import fastinterval

file('README.rst', 'w').write(fastinterval.__doc__)


setup(
    name='fastinterval',
    version=fastinterval.VERSION,
    py_modules=['fastinterval'],
    author='James Casbon',
    author_email='casbon@gmail.com',
    description='Interval class and fasta access',
    long_description=fastinterval.__doc__,
    test_suite='nose.collector',
    requires=['pyfasta', 'bxpython'],
    url='https://github.com/jamescasbon/fastinterval',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
      ],
    keywords='bioinformatics',
)
