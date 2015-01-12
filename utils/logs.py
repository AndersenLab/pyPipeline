from datetime import datetime
from utils import *

class EAV:
    """
    Very simple Entity-Attribute-Value Object
    """

    def __init__(self):
        self.entity = ""
        self.sub_entity = ""
        self.attribute = ""
        self.sub_attribute = ""
        self.value = ""
        self.timestamp = datetime.now()
        self.comment = ""
        self.file = None

    def __repr__(self):
        return "\nEntity:{self.entity}\n\
                Entity:{self.sub_entity}\n\
                Attribute:{self.attribute}\n\
                Sub-Attribute:{self.sub_attribute}\n\
                Value:{self.value}\n\
                timestamp:{self.timestamp}\n".format(**locals())

    def save(self):
        if self.file is None:
            raise Exception("No Log File Set")
        if not file_exists(self.file):
            write_header = True
        else:
            write_header = False
        with(open(self.file, "a")) as f:
            if write_header is True:
                f.write("entity\tsub_entity\tattribute\tsub_attribute\tvalue\tcomment\ttimestamp\n")
            line = '\t'.join(map(str,[self.entity,
                              self.sub_entity,
                              self.attribute,
                              self.sub_attribute,
                              self.value,
                              self.comment,
                              self.timestamp]))
            f.write(line + "\n")


class general_log:
    """
        Simple log file.
    """
    def __init__(self, config):
        self.log = open(config.config_name + ".log", 'a')

    def add(self, msg, analysis_type):
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_msg = "{d}\t{a}\t{l}\n".format(d=now, a=analysis_type, l=msg)
        self.log.write(log_msg)


class command_log:
    """
    The Command log is used to log commands that have been run.
    """
    def __init__(self, config):
        self.log = open(config.config_name + ".commands.log", 'a')

    def add(self, command):
        # Clean up whitespace.
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        command = re.sub("[^\S\r\n]+", " ", command).replace("\n ", "\n").strip()
        log_msg = "{d}\t{comm}\n".format(d=now, comm=command)
        self.log.write(log_msg)
