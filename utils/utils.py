import os
import sys
from subprocess import Popen, PIPE


#======================#
# Set System Specifics #
#======================#
def msg(text, msg_type="\033[1m"):
    """ Reports an error to the user and exits """
    color = {'error': '\033[91m', 'warning': '\033[93m'}
    print("")
    print(color[msg_type] + "pyPipeline %s: " % (msg_type) + '\033[0m' + text)
    print("")
    if msg_type == "error":
        sys.exit(0)


class dotdictify(dict):
    """
        http://stackoverflow.com/questions/3031219
        Creates dictionary accessible via attributes or keys.
    """
    marker = object()

    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise Exception('expected dict')

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, dotdictify):
            value = dotdictify(value)
        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        found = self.get(key, dotdictify.marker)
        if found is dotdictify.marker:
            found = dotdictify()
            dict.__setitem__(self, key, found)
        return found

    __setattr__ = __setitem__
    __getattr__ = __getitem__


def construct_filters(filter_list, soft=True):
    """ Constructs set of piped filters """
    if filter_list:
        filter_command = []
        for k, v in filter_list.items():
            if soft is True:
                filter_command.append("bcftools filter -O u --soft-filter {k} --exclude '{v}'".format(**locals()))
            else:
                filter_command.append("bcftools filter -O u --exclude '{v}'".format(**locals()))
        return '| ' + ' | '.join(filter_command) + ' | bcftools view -O z '
    else:
        return ''


def format_command(command_config):
    """
        Performs standard formatting of commands being run.
    """
    opts = ""
    if command_config is not None:
        for k, v in command_config.items():
            # Use '__' for custom options
            if k.startswith("__"):
                pass
            # Use '_' to designate flags.
            elif k.startswith("_"):
                opts += " %s " % v
            else:
                opts += "%s %s " % (k, v)
        return opts
    else:
        return ""


def file_exists(filename):
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False


def makedir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def is_dir(dir):
    return os.path.isdir(dir)


def remove_file(file):
    try:
        os.remove(file)
    except:
        pass


def which(program):
    """ returns the path to an executable or None if it can't be found
     http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
     """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def version(program):
    """ Attempts to return the version of an executable or None if it can't be found. """
    prog = which(program)
    if program not in ["tabix", "bwa", "rm", "cp", "ls"]:
        try:
            version = Popen([prog, "--version"], stdout=PIPE).communicate()[0]
            version = version.strip().split("\n")[0]
            if version is None:
                version = ""
            return "%-50s\t%s" % (prog, version)
        except:
            return ""
    else:
        # Hand special cases
        return "%-50s\t%s" % (prog, "")


def common_prefix(strings):
    """ Find the longest string that is a prefix of all the strings.
     http://bitbucket.org/ned/cog/src/tip/cogapp/whiteutils.py
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix


def boolify(s):
    """ http://stackoverflow.com/questions/7019283/automatically-type-cast-parameters-in-python """
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")


def set_type(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s


def rreplace(s, old, new, count):
    """
    Replaces last occurance of something
    stackoverflow.com/questions/2556108/
    """
    return (s[::-1].replace(old[::-1], new[::-1], count))[::-1]


def get_column(filename, col_num, delim="\t"):
    """ Returns list of a column from a file """
    column = []
    if file_exists(filename):
        with open(filename, 'r') as r:
            for row in r:
                column.append(row.split(delim)[col_num])
        return column
    else:
        return []
