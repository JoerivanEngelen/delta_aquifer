from setuptools import setup

setup(name='delta_aquifer',
      version='0.1',
      description='Create synthetic delta aquifers based on a few parameters',
      author='Joeri van Engelen',
      author_email='joeri.vanengelen@deltares.nl',
      license='MIT',
      packages=['delta_aquifer'],
      install_requires=[
          'numpy', 'scipy', 'xarray', 'matplotlib'
      ],
      zip_safe=False)