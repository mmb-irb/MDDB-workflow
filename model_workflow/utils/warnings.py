import warnings

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s:%s: %s:\n%s\n' % (filename, lineno, category.__name__, message)

warnings.formatwarning = warning_on_one_line

def warn(msg,category=None):
    warnings.warn(msg,category,stacklevel=2)
    