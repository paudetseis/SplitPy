project('splitpy', 'c',
  version : '0.2.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
)

py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

py.install_sources(
  'splitpy/classes.py',
  'splitpy/calc.py',
  'splitpy/utils.py',
  'splitpy/gui.py',
  'splitpy/__init__.py',
  'splitpy/scripts/split_calc_auto.py',
  'splitpy/scripts/split_calc_manual.py',
  'splitpy/scripts/split_average.py',
  'splitpy/scripts/__init__.py',
  'splitpy/tests/get_meta.py',
  'splitpy/tests/test_args.py',
  'splitpy/tests/test_imports.py',
  'splitpy/tests/test_scripts.py',
  'splitpy/tests/__init__.py',
  pure: false,
  preserve_path: true,
  )
  
  install_data(
    'splitpy/examples/data/MMPY.pkl',
    install_dir: py.get_install_dir(pure: false) / 'splitpy/examples/data',
)
