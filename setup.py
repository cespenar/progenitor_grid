from setuptools import setup

setup(name='progenitor_grid',
    version='0.1.1',
    description='Tool for processing a calculated grid of MESA sdB progenitors.',
    url='https://github.com/cespenar/progenitor_grid',
    author='Jakub Ostrowski',
    author_email='cespenar1@gmail.com',
    license='MIT',
    packages=['progenitor_grid'],
    install_requires=['mesa_reader'])
