#!/usr/bin/env python

from distutils.core import setup

setup(name             = 'oasa',
      version          = '0.14.0b1',
      description      = "OASA is a free cheminformatics library written in Python",
      author           = "Beda Kosata",
      author_email     = "beda@zirael.org",
      maintainer       = "Reinis Danne",
      maintainer_email = "rei4dan@gmail.com"
      url              = "http://bkchem.zirael.org/oasa_en.html",
      license          = "GNU GPL",
      platforms        = ["Unix", "Windows", "hopefully other OSes able to run Python"],

      long_description = "OASA is a free cheminformatics library written in Python",

      packages=['oasa', 'oasa/graph'],
     )

