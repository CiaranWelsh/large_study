from setuptools import setup


#version
MAJOR = 0
MINOR = 0
MICRO = 0



###test version
#MAJOR = 0
#MINOR = 0
#MICRO = 6




__version__ = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

setup(
  name = 'large_study',
  packages = ['large_study'], # this must be the same as the name above
  version = __version__,
  author = 'Ciaran Welsh',
  package_data={'large_study':['*.py']},
  author_email = 'c.welsh2@newcastle.ac.uk',
  url = 'https://github.com/CiaranWelsh/large_study', # use the URL to the github repo
  license='MIT'
)





