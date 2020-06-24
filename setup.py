import os
from setuptools import setup


current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), 'r') as fid:
	LONG_DESCRIPTION = fid.read()

setup(
	name="ImReP",
	version="0.9",
	author="Igor Mandric, Jaqueline Brito, Serghei Mangul",
	description=("ImReP is a method for rapid and accurate profiling of the adaptive immune repertoires from regular RNA-Seq data."),
	long_description=LONG_DESCRIPTION,
	keywords="imrep",
	url="https://github.com/Mangul-Lab-USC/imrep",
	packages=['ImReP'],
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
	install_requires=[
		'argparse',
		'collections-extended',
		'setuptools>=24.2.0',
		'pysam',
		'biopython',
		'intervaltree',
		'jellyfish',
		'numpy',
        'networkx'
	],
	zip_safe=False,
    scripts={'bin/imrep'},
	package_data={'ImReP': ['db/human/*.faa', 'db/mouse/*.faa']},
	classifiers=[
        "Development Status :: 3 - Alpha",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
		"Intended Audience :: Science/Research",
		"Programming Language :: Python :: 3.5",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Natural Language :: English",
		"Operating System :: MacOS :: MacOS X",
		"Operating System :: POSIX :: Linux",
	],
)
