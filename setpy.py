try
    from setuptools import setup
except ImportError:
    from distutils.core import setup

cofig = {
        'description' : 'anuga-cuda project',
        'author' : 'Zhe Weng',
        'url' : 'http://code.google.com/p/anuga-cuda/',
        'download_url' : 'svn checkout http://anuga-cuda.googlecode.com/svn/trunk/ anuga-cuda-read-only',
        'author_email' : 'u5044856@anu.edu.au',
        'version' : '0.1',
        'install_requires' : ['nose'],
        'packages' : ['Boost', 'ANUGA', 'Numpy', 'PyCUDA', 'CUDA'],
        'scripts': [],
        'name' : 'anuga-cuda'
        }

setup(**config)
