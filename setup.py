from setuptools import setup, find_packages

with open("README.rst", "r") as fh:
	long_description = fh.read()

setup(name="pmt_waveform", version=1.0, 
      package_dir={"": "lib"},
      packages=find_packages(), 
      author="Charles Blakemore, Gautam Venugopalan", 
      author_email="chas.blakemore@gmail.com",
      description="PMT Waveform Analysis Library",
      long_description=long_description,
      url="https://github.com/charlesblakemore/pmt_waveform")

