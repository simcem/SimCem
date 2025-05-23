#!/usr/bin/env python3

setup(
    name="simcem",
    version="0.0.2",
    author="Marcus Bannerman",
    author_email="m.bannerman@gmail.com",
    description="SimCem, the thermodynamics/process engineering package",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    url="http://simcem.com",
    classifiers=[
        "Programming Language :: C++",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages('pysrc'),
    package_dir={"":"pysrc"},
    package_data={
        'simcem': ['*.xml'],
    },
    ext_modules=[
        CMakeExtension('simcem.core')
    ],
    install_requires = [
        'cmake',
        'numpy',
        'scipy',
    ],

    test_suite='tests',
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
