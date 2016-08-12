from setuptools import setup

setup(name='ffeatools',
      version='0.9',
      description='FFEA file generation and analysis tools',
      url='http://ffea.bitbucket.com',
      author='FFEA Team',
      author_email='???',
      license='???',
      packages=['ffeatools'],
      install_requires=[
          'numpy',
          'pymol',
          'matplotlib',
          'argparse',
          'math'
      ],
      zip_safe=False)