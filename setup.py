#!/usr/bin/env python
#
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of pytaxonkit (http://github.com/bioforensics/pytaxonkit)
# and is licensed under the BSD license.
#
#
# This Software was prepared for the Department of Homeland Security
# (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
# part of contract HSHQDC-15-C-00064 to manage and operate the National
# Biodefense Analysis and Countermeasures Center (NBACC), a Federally
# Funded Research and Development Center.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -------------------------------------------------------------------------------------------------

from setuptools import setup


with open('README.md', 'r') as infile:
    longdesc = infile.read()

setup(
    name='pytaxonkit',
    description='Python bindings for the TaxonKit library',
    long_description=longdesc,
    long_description_content_type='text/markdown',
    url='https://github.com/bioforensics/pytaxonkit',
    author='Daniel Standage',
    author_email='daniel.standage@nbacc.dhs.gov',
    py_modules=['pytaxonkit'],
    python_requires='>=3.6',
    classifiers=[
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=True,
)
