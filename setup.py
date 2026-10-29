#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re

from setuptools import setup

VERSION = re.findall(
    r"__version__ = \"(.*?)\"", open(os.path.join("desidescsn", "__init__.py")).read()
)[0]


# Where the magic happens:
setup(
    version=VERSION,
)
