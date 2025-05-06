# Copyright (c) 2018      Frederik Beaujean
# Copyright (c) 2017-2025 Danny van Dyk
# Copyright (c) 2021      Philip Lueghausen
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

"""
EOS is a software package that addresses several use cases in the field of
high-energy flavor physics (HEP).

Find documentation online at eos.github.io
Visit us at github.com/eos/eos

"""

from .config import *

# make sure that EOS_HOME points to the location of the wheel supplied data
# if unset.
import os as _os
try:
    if is_wheel:
        if not 'EOS_HOME' in _os.environ:
            _os.environ['EOS_HOME'] = _os.path.normpath(_os.path.join(_os.path.dirname(__file__), '..', '_eos_data/'))
except NameError:
    pass

from _eos import *
from _eos import __version__
from .data import *
from .plot import *
from .datasets import DataSets
from .analysis import Analysis, BestFitPoint
from .analysis_file import AnalysisFile
from .constraint import Constraints
from .ipython import __ipython__
from .observable import Observables
from .parameter import Parameters
from .reference import References
from .signal_pdf import SignalPDF, SignalPDFs
from .tasks import *
from .pyhf_likelihood import PyhfLogLikelihood

import logging
logging.addLevelName(logging.INFO + 1, 'INPROGRESS')
setattr(logging, 'INPROGRESS', logging.INFO + 1)
logging.addLevelName(logging.INFO + 2, 'COMPLETED')
setattr(logging, 'COMPLETED', logging.INFO + 2)
logging.addLevelName(logging.WARNING - 1, 'SUCCESS')
setattr(logging, 'SUCCESS', logging.WARNING - 1)
logger = logging.getLogger('EOS')
logger.setLevel(logging.DEBUG)
logger.propagate = False
# log to stderr by default
stderr_handler = logging.StreamHandler(stream=sys.stderr)
stderr_handler.setLevel(logging.INFO)
logger.addHandler(stderr_handler)

from _eos import _register_log_callback, _set_native_log_level, _NativeLogLevel
_set_native_log_level(_NativeLogLevel.INFO) # default native log level

_MAP_PYTHON_TO_NATIVE_LOG_LEVEL = {
    logging.DEBUG: _NativeLogLevel.DEBUG,
    logging.INFO: _NativeLogLevel.INFO,
    logging.INPROGRESS: _NativeLogLevel.INPROGRESS,
    logging.COMPLETED: _NativeLogLevel.COMPLETED,
    logging.SUCCESS: _NativeLogLevel.SUCCESS,
    logging.WARNING: _NativeLogLevel.WARNING,
    logging.ERROR: _NativeLogLevel.ERROR,
}
def set_log_level(level):
    logger.setLevel(level)
    if level not in _MAP_PYTHON_TO_NATIVE_LOG_LEVEL:
        raise RuntimeError(f'Cannot handle unknown log level: {level}')
    _set_native_log_level(_MAP_PYTHON_TO_NATIVE_LOG_LEVEL[level])

def debug(msg, *args, **kwargs):
    logger.debug(msg, *args, **kwargs)

def info(msg, *args, **kwargs):
    logger.info(msg, *args, **kwargs)

def inprogress(msg, *args, **kwargs):
    logger.log(logging.INPROGRESS, msg, *args, **kwargs)

def completed(msg, *args, **kwargs):
    logger.log(logging.COMPLETED, msg, *args, **kwargs)

def success(msg, *args, **kwargs):
    logger.log(logging.SUCCESS, msg, *args, **kwargs)

def warn(msg, *args, **kwargs):
    logger.warning(msg, *args, **kwargs)

def error(msg, *args, **kwargs):
    logger.error(msg, *args, **kwargs)

def _log_callback(id, level, msg):
    full_msg = f'{id} {msg}'

    if   level == _NativeLogLevel.DEBUG:
        debug(full_msg)
    elif level == _NativeLogLevel.INFO:
        info(full_msg)
    elif level == _NativeLogLevel.INPROGRESS:
        inprogress(full_msg)
    elif level == _NativeLogLevel.COMPLETED:
        completed(full_msg)
    elif level == _NativeLogLevel.SUCCESS:
        success(full_msg)
    elif level == _NativeLogLevel.WARNING:
        warn(full_msg)
    elif level == _NativeLogLevel.ERROR:
        error(full_msg)
    elif level == _NativeLogLevel.SILENT:
        pass
    else:
        raise RuntimeError(f'Cannot handle log level: {level}. Log message: {full_msg}')

_register_log_callback(_log_callback)

import time as _time
import os as _os
def installation_time():
    return _time.ctime(_os.path.getmtime(eos.__file__))

def installation_dir():
    return _os.path.dirname(eos.__file__)

# setup IPython integration
if __ipython__:
    ip = get_ipython()
    html_formatter = ip.display_formatter.formatters['text/html']

    # add the html formatter to the html display hook
    from .ipython import __format_Parameter, __format_KinematicVariable, __format_Kinematics, __format_Options
    from .ipython import __format_Observable, __format_ObservableEntry, __format_GoodnessOfFit, __format_Reference
    html_formatter.for_type(Parameter, __format_Parameter)
    html_formatter.for_type(KinematicVariable, __format_KinematicVariable)
    html_formatter.for_type(Kinematics, __format_Kinematics)
    html_formatter.for_type(Options, __format_Options)
    html_formatter.for_type(Observable, __format_Observable)
    html_formatter.for_type(ObservableEntry, __format_ObservableEntry)
    html_formatter.for_type(GoodnessOfFit, __format_GoodnessOfFit)
    html_formatter.for_type(Reference, __format_Reference)
