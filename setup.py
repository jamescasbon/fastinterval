from setuptools import setup


# file('README.rst', 'w').write(fastinterval.__doc__)


# get the version without an import
VERSION = "Undefined"
DOC = ""
inside_doc = False
for line in open('fastinterval.py'):
    if "'''" in line:
        inside_doc = not inside_doc
    if inside_doc:
        DOC += line.replace("'''", "")

    if (line.startswith('VERSION')):
        exec(line.strip())


setup(
    name='fastinterval',
    py_modules=['fastinterval'],
    author='James Casbon',
    author_email='casbon@gmail.com',
    description='Interval class and fasta access',
    long_description=DOC,
    install_requires=['numpy', 'pyfasta', 'bx-python'],
    url='https://github.com/jamescasbon/fastinterval',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
      ],
    keywords='bioinformatics',
    test_suite='nose.collector',
    version=VERSION
)
