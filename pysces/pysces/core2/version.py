MAJOR = 0
MINOR = 8
MICRO = 0
STATUS = 'rc1'

def current_version():
    return '%s.%s.%s%s' % (MAJOR, MINOR, MICRO, STATUS)

def current_version_tuple():
    return (MAJOR, MINOR, MICRO)

__version__ = current_version() 