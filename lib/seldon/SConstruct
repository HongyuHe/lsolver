import distutils.sysconfig, os
env = Environment(ENV = os.environ,
                  SWIGFLAGS = ['-Wall', '-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")
env.Append(CPPFLAGS = "-DSELDON_DEBUG_LEVEL_4")

conf = Configure(env)
# Link to the appropriate version of Python.
python_version = distutils.sysconfig.get_python_version()
conf.CheckLib("python" + python_version)

if env['PLATFORM'] == 'win32':
	env.Replace(SHLIBSUFFIX = ".pyd")
	env.SharedLibrary('_seldon', ['Seldon.cpp', 'seldon.i'])
else:
	env.SharedLibrary('_seldon.so', ['Seldon.cpp', 'seldon.i'])
