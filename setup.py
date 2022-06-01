from setuptools import setup, find_packages

setup(
    name='MOKIT',
    packages=find_packages(),
    ext_modules= ,
    cmdclass= ,
    license=None,
    description='''Molecular Orbital KIT''',
    author='Jingxiang Zou',
    author_email='njumath@sina.cn',
    url='https://gitlab.com/jxzou/mokit',
    install_requires=[
        "mkl==2019",
        "mkl-include",
        "numpy",
        "gfortran"
        ],
    scripts=["bin/automr"]
    )

