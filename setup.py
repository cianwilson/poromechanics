import setuptools
  
with open("README.md", "r") as fh:
    description = fh.read()
  
setuptools.setup(
    name="poromechanics",
    version="0.0.2",
    author="Matt_McLean",
    author_email="mclean.l.matthew@gmail.com",
    packages=["poromechanics"],
    description="Coupled geomechanics simulation",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/Matt-L-McLean/poromechanics",
    install_requires=[]
)
