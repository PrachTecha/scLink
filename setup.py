from setuptools import setup, find_packages

setup(
    name='scLink',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'anndata',
        'scanpy',
        'scipy',
        'matplotlib',
        'rpy2',
    ],
    author='Prach Techameena',
    author_email='prach.techa@gmail.com',
    description='A package for analyzing link communities in single-cell AnnData objects.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/prachsk/scLink',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)
