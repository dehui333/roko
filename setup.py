#!/usr/bin/env python3
# encoding: utf-8

from distutils.core import setup, Extension


def configuration(parent_package='', top_path=None):
      import numpy
      from numpy.distutils.misc_util import Configuration
      from numpy.distutils.misc_util import get_info

      config = Configuration('',
                             parent_package,
                             top_path)
      config.add_extension('gen',
                           ['gen.cpp', 'generate.cpp', 'models.cpp', 'align_ins.cpp'],
                           extra_objects=['Dependencies/htslib-1.9/libhts.a', '/home/lindehui/miniconda3/envs/py37/lib/libdeflate.so'],
                           extra_compile_args=['-std=c++14'], language='c++',
                           extra_link_args=['-lz', '-lz', '-lm', '-lpthread', '-lcurl', '-lcrypto', '-ledlib'],
                           include_dirs=['Dependencies/htslib-1.9/htslib', 'include'])

      return config

if __name__ == "__main__":
      from numpy.distutils.core import setup

      setup(
        name='gen',
        version='0.0.1',
        configuration=configuration)
