from setuptools import setup, find_packages

setup(name='PyCFMID',
      version='0.0.2',
      description="Python interface to CFM-ID",
      license='MIT',
      author='Hongchao Ji',
      author_email='ji.hongchao@foxmail.com',
      url='https://github.com/hcji/PyCFMID',
	  long_description_content_type="text/markdown",
      packages=find_packages(),
	  install_requires=['requests', 'pubchempy', 'bs4'],
	  include_package_data = True,
	  classifiers=[
	  'Development Status :: 4 - Beta',
	  'Programming Language :: Python :: 3.6',
	  'Programming Language :: Python :: 3.7'
	  ]
     )