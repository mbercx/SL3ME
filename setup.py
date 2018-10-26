from setuptools import setup, find_packages

setup(
    name="SL3ME",
    version="0.1",
    packages=find_packages(exclude=["docs"]),
    install_requires=[
        "pymatgen",
        "numpy",
        "scipy",
        "click"
    ],
    entry_points='''
        [console_scripts]
        slme=SL3ME.cli.cli:main
    '''
)