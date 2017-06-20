#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# @author Veawor Liu
#
# Copy this build script to boost root directory and execute. A reference copy can be downloaded from
# https://dl.bintray.com/boostorg/release/1.64.0/source/
#


import glob
import logging
import os
import os.path
import shutil
import sys
import traceback

from sys import platform as _platform


class SingleLevelFilter(logging.Filter):
    def __init__(self, lvl):
        super().__init__()
        self.lvl = lvl

    def filter(self, record):
        return record.levelno != self.lvl


def create_logger():
    formatter = logging.Formatter('[%(levelname)s][%(filename)s:%(funcName)s():%(lineno)s] %(message)s')
    f1 = SingleLevelFilter(logging.CRITICAL)
    f2 = SingleLevelFilter(logging.ERROR)
    f3 = SingleLevelFilter(logging.WARNING)
    f4 = SingleLevelFilter(logging.INFO)
    f5 = SingleLevelFilter(logging.DEBUG)
    h1 = logging.StreamHandler(sys.stdout)
    h1.addFilter(f1)
    h1.addFilter(f2)
    h1.addFilter(f3)
    h1.setFormatter(formatter)
    h2 = logging.StreamHandler(sys.stderr)
    h2.addFilter(f4)
    h2.addFilter(f5)
    h2.setFormatter(formatter)
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(h1)
    logger.addHandler(h2)
    return logger


log = create_logger()
root_dir = os.path.normpath(os.path.dirname(os.path.realpath(__file__)))
build_dir = os.path.join(root_dir, '_build')


def copy_files(src, dst):
    for file in src:
        if os.path.isfile(file):
            shutil.copy2(file, dst)


def execute(cmd):
    if os.system('{}'.format(cmd)) != 0:
        raise Exception('Failed to execute: {}'.format(cmd))


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def rmdir(path):
    if os.path.exists(path):
        shutil.rmtree(path)


def remove(path):
    if os.path.exists(path):
        os.remove(path)


def build():
    rmdir(build_dir)

    log.info('Building tools...')
    if _platform == "win32":
        execute(os.path.join(root_dir, 'bootstrap'))
    else:
        execute(os.path.join(root_dir, 'bootstrap.sh'))
    b2 = os.path.join(root_dir, 'b2')
    execute('{} {}'.format(b2, os.path.join(root_dir, 'tools', 'bcp')))

    log.info('Exporting headers...')
    mkdir(os.path.join(build_dir, 'include'))
    cmd = os.path.join('dist', 'bin', 'bcp')
    cmd += ' iostreams'
    cmd += ' program_options'
    cmd += ' random'
    cmd += ' serialization'
    cmd += ' {}'.format(os.path.join(build_dir, 'include'))
    execute(cmd)

    if _platform == "win32":
        log.info('Building for win32...')
        cmd = b2
        cmd += ' --with-iostreams --with-program_options --with-random --with-serialization'
        cmd += ' --toolset=msvc-14.1'
        cmd += ' address-model=32 architecture=x86 link=static threading=multi runtime-link=shared'
        cmd += ' --build-type=minimal stage'
        cmd += ' --stagedir={}'.format(os.path.join(root_dir, 'build', 'win32'))
        execute(cmd)
        log.info('Building for x64...')
        cmd = b2
        cmd += ' --with-iostreams --with-program_options --with-random --with-serialization'
        cmd += ' --toolset=msvc-14.1'
        cmd += ' address-model=64 architecture=x86 link=static threading=multi runtime-link=shared'
        cmd += ' --build-type=minimal stage'
        cmd += ' --stagedir={}'.format(os.path.join(root_dir, 'build', 'x64'))
        execute(cmd)
        mkdir(os.path.join(build_dir, 'lib', 'windows', 'win32'))
        mkdir(os.path.join(build_dir, 'lib', 'windows', 'x64'))
        copy_files(glob.glob(os.path.join(root_dir, 'build', 'win32', 'lib', '*.lib')),
                   os.path.join(build_dir, 'lib', 'windows', 'win32'))
        copy_files(glob.glob(os.path.join(root_dir, 'build', 'x64', 'lib', '*.lib')),
                   os.path.join(build_dir, 'lib', 'windows', 'x64'))
    else:
        raise Exception('The {} is not supported.'.format(_platform))

    log.info('Removing intermediate files...')
    rmdir(os.path.join(build_dir, 'include', 'libs'))
    rmdir(os.path.join(build_dir, 'include', 'doc'))
    remove(os.path.join(build_dir, 'include', 'boost.css'))
    remove(os.path.join(build_dir, 'include', 'boost.png'))
    remove(os.path.join(build_dir, 'include', 'Jamroot'))


if __name__ == '__main__':
    try:
        build()
    except Exception as e:
        log.error(e)
        log.error(traceback.format_exc())
