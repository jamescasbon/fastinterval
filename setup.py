from setuptools import setup

setup(
    name='fastinterval',
    py_modules=['fastinterval'],
    author='James Casbon',
    author_email='casbon@gmail.com',
    description='Interval class and fasta access',
    #long_description=fastinterval.__doc__,
    #test_suite='test',
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
