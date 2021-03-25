#!/usr/bin/env python3
# encoding: utf-8


import subprocess
import platform
import pathlib
import os


from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
  def __init__(self, name):
        super().__init__(name=name, sources=[])


class CMakeBuild(build_ext):
  def run(self):
    self.__check_cmake()
    for extension in self.extensions:   
      if isinstance(extesnion, CMakeExtension):
        self.__run_cmake(extension)              

  @staticmethod
  def __check_cmake():
    try:        
      subprocess.check_call(['which', 'cmake'])
      return True
    except subprocess.CalledProcessError:                
      raise RuntimeError('CMake not found')


  def __run_cmake(self, extension):
    # TODO: account for modularity
    # make -H./ extension dependent
    build_dir = pathlib.Path(self.build_temp)
    os.makedirs(build_dir, exist_ok=True)

    cmake_args = [
      '-H./',
      '-B{}'.format(build_dir),
      '-DCMAKE_EXPORT_COMPILE_COMMANDS=1'
    ]

    self.spawn(['cmake'] + cmake_args)
    self.spawn(['cmake', '--build', build_dir])

if __name__ == "__main__":
  if platform.system() == 'Windows':
    raise RuntimeError('Unsupported platform.')
        
  setup(
    name='roko',
    version='0.1.0',
    ext_modules=[CMakeExtension('_roko')],
    cmdclass={'build_ext': CMakeBuild}
  )