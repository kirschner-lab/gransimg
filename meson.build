project(
  'gransimg',
  'c',
  meson_version : '>= 1.3.0',
  version : '0.1',
  default_options : ['warning_level=3'],
  license : 'Apache-2.0',
)

cc = meson.get_compiler('c')

dependencies = [
  cc.find_library('m'),
  dependency('openmp'),
  dependency('libtiff-4', fallback : 'libtiff'),
  dependency('libxml-2.0', fallback : 'libxml2'),
]

sources = [
  'gransimg.c',
]

exe = executable(
  'gransimg',
  sources,
  dependencies : dependencies,
  install : true,
)

test('print-help', exe, args : ['--help'])
