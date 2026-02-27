from setuptools import setup, find_packages

with open('requirements.txt') as f:
  requirements = f.read().splitlines()

setup(
  name = 'rbc2',
  packages = find_packages(include=['rbc2', 'rbc2.*']),
  include_package_data=True,
  install_requires=requirements,
  version = '2026.02.27',
  license='',
  description = 'RetroBioCat 2.0',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = '',
  download_url = '',
  keywords = ['enzyme'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 3'],
)