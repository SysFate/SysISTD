'''
Created on 29 Nov 2022

@author: Francois STÜDER
'''

from setuptools import setup, find_packages

from SysISTD import __version__

setup(
    name='SysISTD',
    version=__version__,
    description='Demultiplexer for SysFate Illumina spatial transcriptomic',
    
    license='GNU Affero General Public License v3.0',

    author='SysFate Lab',
    author_email='sysfatelab@gmail.com',
    url='https://github.com/SysFate/SysISTD',
    
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        "Intended Audience :: Science/Research",
        'License :: OSI Approved :: GNU Affero General Public License v3.0',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Text Processing',
        'Topic :: Text Processing :: General',
    ],
    
    python_requires='>=3.6',
    install_requires=["regex","pandas","biopython"],
    packages=find_packages(),
    py_modules=['SysISTD'],
    
    entry_points={
        'console_scripts': [
            'SysISTD = SysISTD.main:main',
        ]
    },
)


"""with open('README.rst') as file:
    long_description = file.read()

setup(
    
    long_description=long_description,
    long_description_content_type='text/x-rst',
    



    package_dir={'regex': 'regex_3'},
    py_modules=['regex.__init__', 'regex.regex', 'regex._regex_core',
     'regex.test_regex'],
    ext_modules=[Extension('regex._regex', [join('regex_3', '_regex.c'),
      join('regex_3', '_regex_unicode.c')])],
)

with open("README.rst", encoding="ascii") as handle:
    readme_rst = handle.read()

setup(
    long_description=readme_rst,
    project_urls={
        "Documentation": "https://biopython.org/wiki/Documentation",
        "Source": "https://github.com/biopython/biopython/",
        "Tracker": "https://github.com/biopython/biopython/issues",
    },
    cmdclass={"test": test_biopython},
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    include_package_data=True,  # done via MANIFEST.in under setuptools
    install_requires=REQUIRES,
)"""